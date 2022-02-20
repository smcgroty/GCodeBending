#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 10:10:14 2022

@author: stefa
"""

import argparse
import sys
import numpy as np
import math
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import re
from collections import namedtuple

Point2D = namedtuple('Point2D', 'x y')
GCodeLine = namedtuple('GCodeLine', 'x y z e f')


#################   USER INPUT PARAMETERS   #########################

def validateSplineInput():
    class ValidatePoints(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if len(values) < 2:
                raise argparse.ArgumentTypeError(f"{self.dest} requires at least 2 points, starting and ending")
            if not all(i > 0 for i in values):
                raise argparse.ArgumentTypeError(f"{self.dest}: all points must be positive")
            setattr(args, self.dest, values)
    return ValidatePoints

parser = argparse.ArgumentParser(description="Fits a gcode model to a user-defined spline")
parser.add_argument("input_file", help="Input gcode file name")
parser.add_argument("output_file", help="Output file name")
parser.add_argument("--x-spline", "-xs", nargs='+', type=float, required=True, action=validateSplineInput(),
                    help="List of spline x coordinates. First value must be in center of model (usually x-axis midpoint). Minimum 2.")
parser.add_argument("--z-spline", "-zs", nargs='+', type=float, required=True, action=validateSplineInput(),
                    help="List of spline z coordinates. Last value must be equal to or greater than model height. Minimum 2.")
parser.add_argument("--layer-height", "-lh", type=float, default=0.3,
                    help="Layer height of the input gcode. Default 0.3")
parser.add_argument("--max-angle","-ma", type=float, default=30,
                    help="Maximum angle printable with your setup. Default 30 degrees")
parser.add_argument("--disc-length", "-dl", type=float, default=0.01,
                    help="Discretization length for the spline length lookup table. Default 0.01")
parser.add_argument("--spline-angle", "-a", type=float, default=-np.pi/6,
                    help="Final spline angle in radians. Should usually be negative. Default -30 degrees ≈ -0.5236 radians")
try:
    arguments = parser.parse_args()
except argparse.ArgumentTypeError as err:
    print(err)
    sys.exit(1)

INPUT_FILE_NAME = arguments.input_file
OUTPUT_FILE_NAME = arguments.output_file
LAYER_HEIGHT = arguments.layer_height
WARNING_ANGLE = arguments.max_angle

SPLINE_X = arguments.x_spline
SPLINE_Z = arguments.z_spline

SPLINE = CubicSpline(SPLINE_Z, SPLINE_X, bc_type=((1, 0), (1, arguments.spline_angle))) #define spline with BC-conditions

DISCRETIZATION_LENGTH = arguments.disc_length

#################   USER INPUT PARAMETERS END  #########################


SplineLookupTable = [0.0]

nx = np.arange(0,SPLINE_Z[-1],1)

xs = np.arange(0,SPLINE_Z[-1],1)
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.plot(SPLINE_X, SPLINE_Z, 'o', label='data')
ax.plot(SPLINE(xs), xs, label="S")
ax.set_xlim(0, 200)
ax.set_ylim(0, 200)
plt.gca().set_aspect('equal', adjustable='box')
# ax.legend(loc='lower left', ncol=2)
plt.show()


def getNormalPoint(currentPoint: Point2D, derivative: float, distance: float) -> Point2D: #claculates the normal of a point on the spline
    angle = np.arctan(derivative) + math.pi /2
    return Point2D(currentPoint.x + distance * np.cos(angle), currentPoint.y + distance * np.sin(angle))

def parseGCode(currentLine: str) -> GCodeLine: #parse a G-Code line
    thisLine = re.compile('(?i)^[gG][0-3](?:\s+x(?P<x>-?[0-9.]{1,15})|\s+y(?P<y>-?[0-9.]{1,15})|\s+z(?P<z>-?[0-9.]{1,15})|\s+e(?P<e>-?[0-9.]{1,15})|\s+f(?P<f>-?[0-9.]{1,15}))*')
    lineEntries = thisLine.match(currentLine)
    if lineEntries:
        return GCodeLine(lineEntries.group('x'), lineEntries.group('y'), lineEntries.group('z'), lineEntries.group('e'), lineEntries.group('f'))

def writeLine(G, X, Y, Z, F = None, E = None): #write a line to the output file
    outputSting = "G" + str(int(G)) + " X" + str(round(X,5)) + " Y" + str(round(Y,5)) + " Z" + str(round(Z,3))
    if E is not None:
        outputSting = outputSting + " E" + str(round(float(E),5))
    if F is not None:
        outputSting = outputSting + " F" + str(int(float(F)))
    outputFile.write(outputSting + "\n")

"""
# legacy - toooo slow!
def onSplineLength(Zheight) -> float: #calculates a new z height if the spline is followed
    return Zheight #for debugging
    discretizationLength = 0.01 #Steps taken to find the spline height
    currentHeight = 0.00
    currentLength = 0.00
    while currentLength < Zheight:
        currentLength += np.sqrt((SPLINE(currentHeight)-SPLINE(currentHeight+discretizationLength))**2 + discretizationLength**2)
        currentHeight += discretizationLength
    return currentHeight
"""

def onSplineLength(Zheight) -> float: #calculates a new z height if the spline is followed
    for i in range(len(SplineLookupTable)):
        height = SplineLookupTable[i]
        if height >= Zheight:
            return i * DISCRETIZATION_LENGTH
    print("Error! Spline not defined high enough!")

def createSplineLookupTable():
    heightSteps = np.arange(DISCRETIZATION_LENGTH, SPLINE_Z[-1], DISCRETIZATION_LENGTH)
    for i in range(len(heightSteps)):
        height = heightSteps[i]
        SplineLookupTable.append(SplineLookupTable[i] + np.sqrt((SPLINE(height)-SPLINE(height-DISCRETIZATION_LENGTH))**2 + DISCRETIZATION_LENGTH**2))



lastPosition = Point2D(0, 0)
currentZ = 0.0
lastZ = 0.0
currentLayer = 0
relativeMode = False
createSplineLookupTable()

with open(INPUT_FILE_NAME, "r") as gcodeFile, open(OUTPUT_FILE_NAME, "w+") as outputFile:
        for currentLine in gcodeFile:
            if currentLine[0] == ";":   #if NOT a comment
                outputFile.write(currentLine)
                continue
            if currentLine.find("G91 ") != -1:   #filter relative commands
                relativeMode = True
                outputFile.write(currentLine)
                continue
            if currentLine.find("G90 ") != -1:   #set absolute mode
                relativeMode = False
                outputFile.write(currentLine)
                continue
            if relativeMode: #if in relative mode don't do anything
                outputFile.write(currentLine)
                # continue
            currentLineCommands = parseGCode(currentLine)
            if currentLineCommands is not None: #if current comannd is a valid gcode
                if currentLineCommands.z is not None: #if there is a z height in the command
                    currentZ = float(currentLineCommands.z)

                if currentLineCommands.x is None or currentLineCommands.y is None: #if command does not contain x and y movement it#s probably not a print move
                    if currentLineCommands.z is not None: #if there is only z movement (e.g. z-hop)
                        outputFile.write("G91\nG1 Z" + str(currentZ-lastZ))
                        if currentLineCommands.f is not None:
                            outputFile.write(" F" + str(currentLineCommands.f))
                        outputFile.write("\nG90\n")
                        lastZ = currentZ
                        continue
                    outputFile.write(currentLine)
                    continue
                currentPosition = Point2D(float(currentLineCommands.x), float(currentLineCommands.y))
                midpointX = lastPosition.x + (currentPosition.x - lastPosition.x) / 2  #look for midpoint

                distToSpline = midpointX - SPLINE_X[0]

                #Correct the z-height if the spline gets followed
                correctedZHeight = onSplineLength(currentZ)

                angleSplineThisLayer = np.arctan(SPLINE(correctedZHeight, 1)) #inclination angle this layer

                angleLastLayer = np.arctan(SPLINE(correctedZHeight - LAYER_HEIGHT, 1)) # inclination angle previous layer

                heightDifference = np.sin(angleSplineThisLayer - angleLastLayer) * distToSpline * -1 # layer height difference

                transformedGCode = getNormalPoint(Point2D(correctedZHeight, SPLINE(correctedZHeight)), SPLINE(correctedZHeight, 1), currentPosition.x - SPLINE_X[0])

                #Check if a move is below Z = 0
                if float(transformedGCode.x) <= 0.0:
                    print("Warning! Movement below build platform. Check your spline!")

                #Detect unplausible moves
                if transformedGCode.x < 0 or np.abs(transformedGCode.x - currentZ) > 50:
                    print("Warning! Possibly unplausible move detected on height " + str(currentZ) + " mm!")
                    outputFile.write(currentLine)
                    continue
                #Check for self intersection
                if (LAYER_HEIGHT + heightDifference) < 0:
                    print("ERROR! Self intersection on height " + str(currentZ) + " mm! Check your spline!")

                #Check the angle of the printed layer and warn if it's above the machine limit
                if angleSplineThisLayer > (WARNING_ANGLE * np.pi / 180.):
                    print("Warning! Spline angle is", (angleSplineThisLayer * 180. / np.pi), "at height  ", str(currentZ), " mm! Check your spline!")

                if currentLineCommands.e is not None: #if this is a line with extrusion
                    """if float(currentLineCommands.e) < 0.0:
                        print("Retraction")"""
                    extrusionAmount = float(currentLineCommands.e) * ((LAYER_HEIGHT + heightDifference)/LAYER_HEIGHT)
                    #outputFile.write(";was" + currentLineCommands.e + " is" + str(extrusionAmount) + " diff" + str(int(((LAYER_HEIGHT + heightDifference)/LAYER_HEIGHT)*100)) + "\n")
                else:
                    extrusionAmount = None
                writeLine(1,transformedGCode.y, currentPosition.y, transformedGCode.x, None, extrusionAmount)
                lastPosition = currentPosition
                lastZ = currentZ
            else:
                outputFile.write(currentLine)
print("GCode bending finished!")
