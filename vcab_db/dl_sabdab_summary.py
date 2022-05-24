#!/usr/bin/env python
# coding: utf-8

import urllib.request, urllib.error, urllib.parse
from bs4 import BeautifulSoup
import wget
from datetime import date

# automatically download Sabdab table of Human Ab structures with constant regions

url = 'http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/search/?ABtype=All&method=All&species=HOMO+SAPIENS&resolution=&rfactor=&antigen=All&ltype=All&constantregion=True&affinity=All&isin_covabdab=All&isin_therasabdab=All&chothiapos=&restype=ALA&field_0=Antigens&keyword_0='
htmlout = urllib.request.urlopen(url).read()

# find the Downloads section and retrieve the url for the summary tsv file
soup = BeautifulSoup(htmlout, 'lxml')
dl = soup.find('div', id='downloads')
download_url = 'http://opig.stats.ox.ac.uk' + dl.findAll('li')[0].findAll('a')[0]['href']

# download the file and save as '<date>_sabdab.tsv'
destination = date.today().strftime("%Y%m%d") + '_sabdab.tsv'
wget.download(download_url, destination)
