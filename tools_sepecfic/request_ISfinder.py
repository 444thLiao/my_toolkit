from bs4 import BeautifulSoup
import re
import pandas as pd

soup = BeautifulSoup(open("/home/liaoth/data2/project/shenzhen_actineto/ISfinder/43458.html").read(), 'html.parser')

body = soup.findAll(name="article")
body = body[0]

all_childrens = list(body.children)[4:]

dom_dict = {}
name = ''
for dom in all_childrens:
    if dom.name == 'pre' and 'spades=' in dom.text:
        name = re.findall('spades=(.*)',dom.text)
        if name:
            name = name[0]
    elif dom.name == 'table':
        empty_df = pd.DataFrame()
        tbody = dom.findNext(name='tbody')
        # tbody_childrens = tbody.children
        for row in tbody.children:
            if row.name:
                empty_df = empty_df.append(pd.DataFrame([[_.text for _ in row.children]],index=[0]))
            else:
                dom_dict[name] = empty_df
