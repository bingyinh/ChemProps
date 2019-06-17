#!/usr/bin/python

# 06/04/2019 Bingyin Hu

from pymongo import MongoClient
import logging

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

    # the main search function for polymer infos
    # will call five sub-methods for mapping, wf stands for weighting factor
    # 1) apple to apple comparison for polymer names (in _stdname, _abbreviations, _synonyms), wf 3
    # 2) apple to apple comparison for abbreviations (in _abbreviations), wf 2
    # 3) relaxed bag-of-word comparison for tradenames (in _tradenames), wf 1
    # 4) bag-of-character comparison for polymer names (in _stdname, _synonyms), wf 3
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
        for cand in self.cp.polymer.find({'_abbreviations': {'$regex': rptname, '$options': 'i'}})
            if cand['_id'] not in candidates:
                candidates[cand['_id']] = {'data': cand, 'wf': 0}
            candidates[cand['_id']]['wf'] += 3
        # query for '_synonyms' array
        for cand in self.cp.polymer.find({'_synonyms': {'$regex': rptname, '$options': 'i'}})
            if cand['_id'] not in candidates:
                candidates[cand['_id']] = {'data': cand, 'wf': 0}
            candidates[cand['_id']]['wf'] += 3
        # 2) apple to apple comparison for abbreviations (in _abbreviations), wf 2
        if 'Abbreviation' in keywords:
            rptabbr = keywords['Abbreviation']
            # query for '_abbreviations' array
            for cand in self.cp.polymer.find({'_abbreviations': {'$regex': rptabbr, '$options': 'i'}})
                if cand['_id'] not in candidates:
                    candidates[cand['_id']] = {'data': cand, 'wf': 0}
                candidates[cand['_id']]['wf'] += 2
        # 3) relaxed bag-of-word comparison for tradenames (in _tradenames), wf 1
        if 'TradeName' in keywords:
            rpttrad = keywords['TradeName']
            # query for '_tradenames' array
            for cand in self.cp.polymer.find({'_tradenames': {'$regex': '/%s/' %(rpttrad), '$options': 'i'}})
        # 4) bag-of-character comparison for polymer names (in _stdname, _synonyms), wf 3
        # 5) relaxed bag-of-word comparison for polymer names (in _stdname, _synonyms), wf 2


    def containAllWords(self, query, field, collection):
        import re
        import string
        pattern = re.compile('[^a-zA-Z]', re.UNICODE)
        pattern.sub(' ', query) 
        words = query.split()
        ids = dict()
        for word in words:
            for result in collection.find({field: {'$regex': word, '$options':'i'}}):
                if result['_id'] not in ids:
                    ids[result['_id']] = {'data': result, 'freq': 0}
                ids[result['_id']]['freq'] += 1
        # check if the length of words equals any of the freq in ids
