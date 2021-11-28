import json
from json import JSONEncoder
import numpy


class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.ndarray):
            return obj.tolist()
        return JSONEncoder.default(self, obj)


def openJsonFile():
    file = open("C:/Users/Daniel/Desktop/inputJson.json", "r")
    return file


def serializeJson(file):
    jsonFile = json.load(file)
    return jsonFile


def saveToJson(dty):
    rawDTY = dty
    dictionary = {
        "plotType": "both",
        "derivativeNumber": 0,
        "scaleType": "exponential",
        "plotRadius": 5,
        "atoms": [
            {
                "atomicNumber": 6,
                "shellOccupation": "2.4",
                "dtyRAW": rawDTY
            }
        ]
    }

    with open("C:/Users/Daniel/Desktop/outputJson.json", "w") as outfile:
        json.dump(dictionary, outfile)


def saveJsonTest(dty):
    encodedNumpyData = json.dumps(dty, cls=NumpyArrayEncoder)  # use dump() to write array into file
    return encodedNumpyData
