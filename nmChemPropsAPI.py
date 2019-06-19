#!/usr/bin/python

# 06/04/2019 Bingyin Hu

from pymongo import MongoClient
import logging
import re # for query reformat
import string # for query reformat

class nmChemPropsAPI():
    def __init__(self, nmid):
        # input the NanoMine ID, example: L999_Someone_2020_S2
        self.nmid = nmid
        # load logger
        logging.basicConfig(filename='nmChemPropsInteract.log',
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            level = logging.INFO
                           )
        self.loadMGconfig()
        # mongo init
        self.client = MongoClient('mongodb://%s:%s@localhost:27017/tracking?authSource=admin'
                                  %(self.env['NM_MONGO_USER'],
                                    self.env['NM_MONGO_PWD']
                                   )
                                 )
        # access to DBs
        self.cp = self.client.ChemProps
        self.uk = self.client.unknowns
        
    # load mongo configurations
    def loadMGconfig(self):
        self.env = dict()
        # read mongo.config for configurations
        with open("mongo.config", "r") as f:
            confs = f.read().split('\n')
        for i in range(len(confs)):
            kv = confs[i]
            k = kv.split(':')[0].strip()
            v = kv.split(':')[1].strip()
            self.env[k] = v

    # a helper method to generate a bocStr for a string based on the occurrence
    # of the chars. Example: (a,b,c,...,y,z,0,1,2,3,4,5,6,7,8,9) to 101214...01
    def bagOfChar(self, myStr):
        bag = []
        for alphabet in string.ascii_lowercase:
            bag.append(str(myStr.lower().count(alphabet)))
        for number in string.digits:
            bag.append(str(myStr.lower().count(number)))
        return ''.join(bag)
    # the main search function for polymer infos
    # will call five sub-methods for mapping, wf stands for weighting factor
    # 1) apple to apple comparison for polymer names (in _stdname, _abbreviations, _synonyms), wf 3
    # 2) apple to apple comparison for abbreviations (in _abbreviations), wf 2
    # 3) relaxed bag-of-word comparison for tradenames (in _tradenames), wf 1
    # 4) bag-of-character comparison for polymer names (in _boc), wf 2+1
    # 5) relaxed bag-of-word comparison for polymer names (in _stdname, _synonyms), wf 2
    # input format:
    #   {'ChemicalName': 'Poly(styrene)', 'Abbreviation': 'PS', 'TradeName': 'Dylite', 'uSMILE': ''}
    #   NanoMine schema guarantees 'ChemicalName' has minimum occurrence of 1
    #   'Abbreviation', 'TradeName', and 'uSMILE' are not required, since users might leave them blank
    # output format:
    # if there is a match:
    #   {'StandardName': _stdname, 'uSMILE': _id, 'density': _density}
    #   for multiple matches, return the one with highest cummulated wf
    #   before exit, examine again whether the reported 'Abbreviation' and 'TradeName' are recorded in ChemProps, log them and manually check, if confirmed to be correct, add them to the google spreadsheet
    # if there is not a match:
    #   insert _inputname, _inputabbr, _inputsmiles, _nmid[] to unknowns.polymer
    def searchPolymers(self, keywords):
        candidates = dict() # use '_id' as keys
        # 1) apple to apple comparison for polymer names (in _stdname, _abbreviations, _synonyms), wf 3
        rptname = keywords['ChemicalName']
        # query for '_stdname' with rptname
        for cand in self.cp.polymer.find({'_stdname': {'$regex': rptname, '$options': 'i'}}):
            if cand['_id'] not in candidates:
                candidates[cand['_id']] = {'data': cand, 'wf': 0}
            candidates[cand['_id']]['wf'] += 3
        # query for '_abbreviations' array
        for cand in self.cp.polymer.find({'_abbreviations': {'$regex': rptname, '$options': 'i'}}):
            if cand['_id'] not in candidates:
                candidates[cand['_id']] = {'data': cand, 'wf': 0}
            candidates[cand['_id']]['wf'] += 3
        # query for '_synonyms' array
        for cand in self.cp.polymer.find({'_synonyms': {'$regex': rptname, '$options': 'i'}}):
            if cand['_id'] not in candidates:
                candidates[cand['_id']] = {'data': cand, 'wf': 0}
            candidates[cand['_id']]['wf'] += 3
        # 2) apple to apple comparison for abbreviations (in _abbreviations), wf 2
        if 'Abbreviation' in keywords:
            rptabbr = keywords['Abbreviation']
            # query for '_abbreviations' array
            for cand in self.cp.polymer.find({'_abbreviations': {'$regex': rptabbr, '$options': 'i'}}):
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0}
                candidates[cand['_id']]['wf'] += 2
        # 3) relaxed bag-of-word comparison for tradenames (in _tradenames), wf 1
        if 'TradeName' in keywords:
            rpttrad = keywords['TradeName']
            # query for '_tradenames' array
            relBOW = self.containAllWords(rpttrad, '_tradenames', self.cp.polymer)
            for cand in relBOW:
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0}
                candidates[cand['_id']]['wf'] += 1
            # query for bag-of-character for only alphabets (use $regex '^0100020...')
            tradnameBOC = self.bagOfChar(rpttrad)
            tradnameBOCalph = tradnameBOC[:-10] # remove number's index
            for cand in self.cp.polymer.find({'_boc': {'$regex': '^%s' %(tradnameBOCalph)}}):
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0}
                candidates[cand['_id']]['wf'] += 1
        # 4) bag-of-character comparison for polymer names (in _boc), wf 2
        rptnameBOC = self.bagOfChar(rptname)
        for cand in self.cp.polymer.find({'_boc': {'$regex': '%s' %(rptnameBOC)}}):
            if cand['_id'] not in candidates:
                candidates[cand['_id']] = {'data': cand, 'wf': 0}
            candidates[cand['_id']]['wf'] += 2
        # only alphabets version, wf 1
        rptnameBOCalph = rptnameBOC[:-10] # remove number's index
        for cand in self.cp.polymer.find({'_boc': {'$regex': '^%s' %(rptnameBOCalph)}}):
            if cand['_id'] not in candidates:
                candidates[cand['_id']] = {'data': cand, 'wf': 0}
            candidates[cand['_id']]['wf'] += 1
        # 5) relaxed bag-of-word comparison for polymer names (in _stdname, _synonyms), wf 2
        # query for '_stdname' array
        relBOWstd = self.containAllWords(rptname, '_stdname', self.cp.polymer)
        for cand in relBOWstd:
            if cand['_id'] not in candidates:
                candidates[cand['_id']] = {'data': cand, 'wf': 0}
            candidates[cand['_id']]['wf'] += 2
        relBOWsyn = self.containAllWords(rptname, '_synonyms', self.cp.polymer)
        for cand in relBOWsyn:
            if cand['_id'] not in candidates:
                candidates[cand['_id']] = {'data': cand, 'wf': 0}
            candidates[cand['_id']]['wf'] += 2
        # end of the query part
        # if there is not a match:
        #   insert _inputname, _inputabbr, _inputsmiles, _nmid[] to unknowns.polymer
        if len(candidates) == 0:
            
    # if there is a match:
    #   {'StandardName': _stdname, 'uSMILE': _id, 'density': _density}
    #   for multiple matches, return the one with highest cummulated wf
    #   before exit, examine again whether the reported 'Abbreviation' and 'TradeName' are recorded in ChemProps, log them and manually check, if confirmed to be correct, add them to the google spreadsheet

        # return the candidate with the highest wf
        
    # this function querys collection.field array for the collections that
    # contain all of the alphabetic words in the query
    # input:
    #   query ------- a string for query, e.g. "Poly[(polytetrahydrofuran 1000)-alt-(4,4'-diphenylmethane diisocynate)]"
    #   field ------- the field to conduct search, e.g. "_tradenames"
    #   collection -- the collection that contains the field, e.g. "self.cp.polymer"
    # output:
    #   candidates -- a list of collection dicts that meet the criteria
    def containAllWords(self, query, field, collection):
        pattern = re.compile('[^a-zA-Z]', re.UNICODE)
        query = pattern.sub(' ', query)
        words = query.split()
        print words
        ids = dict()
        for word in words:
            for result in collection.find({field: {'$regex': word, '$options':'i'}}):
                if result['_id'] not in ids:
                    ids[result['_id']] = {'data': result, 'freq': 0}
                ids[result['_id']]['freq'] += 1
        # check if the length of words equals any of the freq in ids
        output = []
        nWords = len(words)
        for cand in ids:
            print cand
            if ids[cand]['freq'] == nWords:
                output.append(ids[cand]['data'])
        return output