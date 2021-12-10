import json


def openJsonFile(path):
    file = open(path, "r")
    return file


def serializeJson(file):
    jsonFile = json.load(file)
    return jsonFile


def saveToJson(atomList):

    dictionary = {
        "plotType": "both",
        "derivativeNumber": 0,
        "scaleType": "exponential",
        "plotRadius": 5,
        "atoms": atomList
    }

    with open("C:/Users/Daniel/Desktop/outputJson.json", "w") as outfile:
        json.dump(dictionary, outfile)
