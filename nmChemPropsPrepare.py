#!/usr/bin/python

import requests
import xlrd
from pymongo import MongoClient
import os

class nmChemPropsPrepare():
    def __init__(self):
        self.loadGSconfig()
        # downloadGS
        self.downloadGS()
        # if we changed something in mongo, gsUpdate list will record these changes
        # admins need to manually address these changes in the google spreadsheet
        self.gsUpdate = []
        # prepare filler and polymer data
        self.filler = dict()
        self.polymer = dict()
        self.prepFiller()
        self.prepPolymer()
        # mongo init
        self.client = MongoClient('mongodb://%s:%s@localhost:27017'
                                  %(os.environ['NM_MONGO_USER'],
                                    os.environ['NM_MONGO_PWD']
                                   )
                                 )

    # load google spreadsheet configurations
    def loadGSconfig(self):
        # read gs.config for configurations
        with open("gs.config", "r") as f:
            confs = f.read().split('\n')
        self.url_format = confs[0]
        self.key = confs[1]
        self.format = "xlsx"
        # gids
        self.gids = dict()
        for i in range(2,len(confs)):
            kv = confs[i]
            k = kv.split(':')[0].strip()
            v = kv.split(':')[1].strip()
            self.gids[k] = v
        
    # download google spreadsheets
    def downloadGS(self):
        for fname in self.gids:
            resp = requests.get(self.url_format %(self.key,
                                                  self.format,
                                                  self.gids[fname]
                                                 )
                               )
            with open(fname + ".xlsx", "wb") as f:
                f.write(resp.content)
                print("%s sheet is downloaded as %s.xlsx" %(fname, fname))

    # prepare ChemProps.polymer data
    def prepPolymer(self):
        xlfile = xlrd.open_workbook("matrixRaw.xlsx")
        sheet = xlfile.sheets()[0] # one sheet per xlsx file
        header = sheet.row_values(0) # SMILES;uSMILES;std_name;density(g/cm3);density_std_err(g/cm3);abbreviations;synonyms;tradenames
        # create a map for headers
        hmap = {}
        for i in range(len(header)):
            hmap[header[i]] = i
        # loop
        for row in range(1, sheet.nrows):
            rowdata = sheet.row_values(row)
            # skip the unfilled items
            if (len(rowdata[hmap['abbreviations']]) == 0 and
                len(rowdata[hmap['synonyms']]) == 0 and
                len(rowdata[hmap['tradenames']]) == 0):
                continue
            # otherwise save the data to self.polymer
            if rowdata[hmap['uSMILES']] not in self.polymer:
                self.polymer[rowdata[hmap['uSMILES']]] = {
                    "_id": rowdata[hmap['uSMILES']],
                    "_stdname": rowdata[hmap['std_name']],
                    "_abbreviations": [],
                    "_synonyms": [],
                    "_tradenames": [],
                    "_density": rowdata[hmap['density(g/cm3)']]
                }
            # abbreviations
            if len(rowdata[hmap['abbreviations']]) > 0:
                self.polymer[rowdata[hmap['uSMILES']]]['_abbreviations'] = self.striplist(rowdata[hmap['abbreviations']].split(';'))
            # synonyms
            if len(rowdata[hmap['synonyms']]) > 0:
                self.polymer[rowdata[hmap['uSMILES']]]['_synonyms'] = self.striplist(rowdata[hmap['synonyms']].split(';'))
            # tradenames
            if len(rowdata[hmap['tradenames']]) > 0:
                self.polymer[rowdata[hmap['uSMILES']]]['_tradenames'] = self.striplist(rowdata[hmap['tradenames']].split(';'))


    # prepare ChemProps.filler data
    def prepFiller(self):
        xlfile = xlrd.open_workbook("fillerRaw.xlsx")
        sheet = xlfile.sheets()[0] # one sheet per xlsx file
        header = sheet.row_values(0) # nm_entry/std_name/density_g_cm3
        # create a map for headers
        hmap = {}
        for i in range(len(header)):
            hmap[header[i]] = i
        # loop
        for row in range(1, sheet.nrows):
            rowdata = sheet.row_values(row)
            if rowdata[hmap['std_name']] not in self.filler:
                self.filler[rowdata[hmap['std_name']]] = {"_id":rowdata[hmap['std_name']], "_density": rowdata[hmap['density_g_cm3']], "_alias":[]}
            self.filler[rowdata[hmap['std_name']]]['_alias'].append(rowdata[hmap['nm_entry']])

    # update MongoDB
    def updateMongoDB(self):
        dbnames = client.list_database_names() # check if db exists
        initPolymer = False # a flag inidicating whether this is the first time creating the ChemProps.polymer collection
        initFiller = False # a flag inidicating whether this is the first time creating the ChemProps.filler collection
        if u'ChemProps' not in dbnames:
            initPolymer = True
            initFiller = True
        cp = client.ChemProps
        clctnames = cp.list_collection_names() # check if collection exists
        # if ChemProps exists
        if not initPolymer and 'polymer' not in clctnames:
            initPolymer = True
        if not initPolymer and 'filler' not in clctnames:
            initFiller = True
        ## first creation cases (polymer)
        if initPolymer:
            pol = cp.polymer
            posted_polymer = pol.insert_many(self.polymer.values())
        ## update cases (polymer)
        else:
            # loop through the items in the self.polymer dict, see if everything matches
            for uSMILES in self.polymer:
                gsData = self.polymer['uSMILES'] # google spreadsheet data
                mgData = cp.polymer.find({"_id": uSMILES})[0] # mongo data, find by _id
                # continue if there is no difference between gsData and mgData
                if gsData == mgData:
                    continue
                # otherwise, find the difference
                # if gsData is a superset of mgData, update mgData
                # if mgData is a superset of gsData, record the difference in self.gsUpdate
                result = self.compareDict(d1 = gsData,
                                          d1name = 'google sheet',
                                          d2 = mgData,
                                          d2name = 'mongo',
                                          imtbKeys = {'_id',
                                                      '_stdname',
                                                      '_density'
                                                     }
                                         )
        ## first creation cases (filler)
        if initFiller:
            fil = cp.filler
            posted_filler = fil.insert_many(self.filler.values())


    # remove leading and trailing white spaces
    def striplist(self, mylist):
        for i in range(len(mylist)):
            mylist[i] = mylist[i].strip()
        return mylist

    # compare two dicts, need to specify the keys to the immutable objects,
    # the function returns a result dict that indicates objects that do not
    # exist in the current dict but exist in the other dict for each dict.
    # DO NOT SUPPORT NESTED DICTS
    # example:
    # d1 = {'k1': [1,2], 'k2': 'new'}
    # d2 = {'k1': [1,3], 'k2': 'old'}
    # result = {'d1': [('add', 'k1', 3)], 'd2': [('add', 'k1', 2), ('modify', 'k2', 'new')]}
    def compareDict(self, d1, d1name, d2, d2name, imtbKeys):
        result = {d1name: [], d2name: []} # init output dict
        # prepPolymer guarantees d1 and d2 will have the same keys set even if
        # some keys will have empty string or list
        allKeys = set(d1.keys())
        for key in allKeys:
            # immutables always trust d1 has the latest version
            if key in imtbKeys:
                if d1[imtbKey] != d2[imtbKey]:
                    result[d2name].append(('modify', imtbKey, d1[imtbKey]))
            else:
                # use set.difference() function to get the result
                
        # non immutables

if __name__ == '__main__':
    nm = nmChemPropsPrepare()
    nm.downloadGS()