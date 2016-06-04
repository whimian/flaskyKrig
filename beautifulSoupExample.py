from bs4 import BeautifulSoup

data = ""
with open("static\images\krigingpurple.svg", 'r') as file:
    for line in file.readlines():
        data += line.rstrip('\n')

soup = BeautifulSoup(data, "html.parser")

print(str(soup.svg))
