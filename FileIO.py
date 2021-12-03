import json



def openJsonFile(path):
    file = open(path, "r")
    return file


def serializeJson(file):
    jsonFile = json.load(file)
    return jsonFile


def saveToJson(dty):
    print("data type being fed into FileIO.saveToJson() is " + str(type(dty)))
    dtyList = dty.tolist()
    print("data type being dumped into json file is " + str(type(dtyList)))
    dictionary = {
        "plotType": "both",
        "derivativeNumber": 0,
        "scaleType": "exponential",
        "plotRadius": 5,
        "atoms": [
            {
                "atomicNumber": 6,
                "shellOccupation": "2.4",
                "dtyRAW": dtyList
            }
        ]
    }

    with open("C:/Users/Daniel/Desktop/outputJson.json", "w") as outfile:
        json.dump(dictionary, outfile)
