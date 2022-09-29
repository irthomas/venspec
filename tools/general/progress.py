# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 13:15:28 2022

@author: iant

https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters
"""


def progress(iterable, length=50):
    """display a progress bar in the console"""

    total = len(iterable)

    def printProgressBar (iteration):
        percent = ("{0:.1f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = "*" * filledLength + '-' * (length - filledLength)
        print(f'\r |{bar}| {percent}% ', end = "\r")

    printProgressBar(0)

    for i, item in enumerate(iterable):
        yield item
        printProgressBar(i + 1)
    print()
    
