import urllib.request
import os

print('Beginning file download with urllib2...')

url = 'https://dnamod.hoffmanlab.org/DNAmod.sqlite'
path = os.getcwd()+'/DNAmod.sqlite'
print(str(path))
urllib.request.urlretrieve(url, path)  