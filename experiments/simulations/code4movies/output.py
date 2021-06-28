import os

def createdir(dirname):

    if not os.path.exists(dirname):
        os.makedirs(dirname)

    return()
