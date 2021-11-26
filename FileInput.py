import json


def openJsonFile():
    file = open("C:/Users/Daniel/Desktop/inputJson.json", "r")
    return file


def serializeJson(file):
    jsonFile = json.load(file)
    # print(jsonFile["atoms"])
    return jsonFile
