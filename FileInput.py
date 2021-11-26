import json


def readInputFile():
    json.load()


def openTextFile():
    file = open("C:/Users/Daniel/Desktop/text.txt", "r")
    return file


def openJsonFile():
    file = open("C:/Users/Daniel/Desktop/inputJson.json", "r")
    return file


def readFile(file):
    data = file.read()
    print(data)


def serializeJson(file):
    jsonFile = json.load(file)
    print(jsonFile)
