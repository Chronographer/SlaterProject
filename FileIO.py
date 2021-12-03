import json
from json import JSONEncoder
import numpy


class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.ndarray):
            return obj.tolist()
        return JSONEncoder.default(self, obj)


def openJsonFile(path):
    file = open(path, "r")
    return file


def serializeJson(file):
    jsonFile = json.load(file)
    return jsonFile


def saveToJson(dty):
    dictionary = {
        "plotType": "both",
        "derivativeNumber": 0,
        "scaleType": "exponential",
        "plotRadius": 5,
        "atoms": [
            {
                "atomicNumber": 6,
                "shellOccupation": "2.4",
                "dtyRAW": dty
            }
        ]
    }

    with open("C:/Users/Daniel/Desktop/outputJson.json", "w") as outfile:
        json.dump(dictionary, outfile, cls=NumpyArrayEncoder)


def saveJsonTest(dty):
    print("data type being fed into FileIO.saveJsonTest() is " + str(type(dty)))
    encodedNumpyData = json.dumps(dty, cls=NumpyArrayEncoder)  # use dump() to write array into file
    return encodedNumpyData
