from Tkinter import *
import tkFileDialog as filedialog
from shapely.geometry import Polygon, mapping, MultiPolygon, shape, Point, LineString
import shapely.wkt, geojson, json
import numpy as np
import pysal as ps
import pandas
import fiona, os, csv, math
import geopandas as gpd
from geopy.distance import vincenty
from shapely.ops import cascaded_union, unary_union, transform
import utm
import random
import webbrowser
import operator
from geopandas.tools import sjoin
import folium
from folium.plugins import MarkerCluster
from folium import IFrame
import unicodedata
from functools import partial
import pyproj
from copy import deepcopy

#A function to transform things from UTM to lat / long
#as folium only works with lat / lon. Specifically, this
#code only applies to our area of Montana. A better alternative
#will need to be found later on if we post this online for whatever
#region a person wants to use it in.
project = partial(
    pyproj.transform,
    pyproj.Proj(init='epsg:32612'),
    pyproj.Proj(init='epsg:4326'))

class file_names:
    fieldFileName = ""
    yieldFileName = ""
    proteinFileName = ""
    appliedFileName = ""
    choice = ""
    complete = False
    satisfied = False
    def __init__(self, string):
        self.val = string

class field_vals:
    numYieldBins = -1
    numProteinBins = -1
    nitrogenList = []
    nitrogenExpectations = {}
    nitrogenIndex = {}
    yieldBounds = []
    proteinBounds = []
    gridWidth = -1
    gridHeight = -1
    bufferGiven = 0
    xmin = -1
    xmax = -1
    ymin = -1
    ymax = -1
    yieldMin = -1
    yieldMax = -1
    proteinMin = -1
    proteinMax = -1
    def __init__(self, string):
        self.val = string

class GridCell:
    #change to nitrogen later on
    protein = -1
    cropYield = -1
    yieldBin = -1
    proteinBin = -1
    nitrogen = -1
    entryNum = -1
    isSorted = 0
    def __init__(self, bl_x, bl_y, ur_x, ur_y):
        self.bl_x = bl_x
        self.bl_y = bl_y
        self.ur_x = ur_x
        self.ur_y = ur_y
        self.trueBounds = Polygon([(bl_x, bl_y), (ur_x, bl_y), (ur_x, ur_y), (bl_x, ur_y)])
        self.rectangle = Polygon([(bl_x + .5, bl_y + .5), (ur_x - .5, bl_y + .5),
        (ur_x - .5, ur_y -.5), (bl_x + .5, ur_y - .5)])

#all of the below functions grab filenames and other
#values from the GUI
def field_file(file_names, input, entryBox):
    file_names.fieldFileName = input
    entryBox.insert(END, input)

def yield_file(file_names, input, entryBox):
    file_names.yieldFileName = input
    entryBox.insert(END, input)

def protein_file(file_names, input, entryBox):
    file_names.proteinFileName = input
    entryBox.insert(END, input)

def applied_file(file_names, input, entryBox):
    file_names.appliedFileName = input
    entryBox.insert(END, input)

def option(file_names, optionVar):
    file_names.choice = optionVar.get()

def yieldBinSelect(field_vals, yieldBinOptionVar):
    field_vals.numYieldBins = yieldBinOptionVar.get()

def proteinBinSelect(field_vals, proteinBinOptionVar):
    field_vals.numProteinBins = proteinBinOptionVar.get()

def nitrogen_list(field_vals, entryBox, input):
    fileName = input
    entryBox.insert(END, input)
    f = open(fileName,'r')
    stringNum = f.read()
    nitrogenList = map(int, stringNum.split())
    print nitrogenList
    field_vals.nitrogenList = nitrogenList
    f.close()

def height_val(field_vals, heightEntry):
    stringNum = heightEntry.get()
    value = int(stringNum)
    field_vals.gridHeight = value

def width_val(field_vals, widthEntry):
    stringNum = widthEntry.get()
    value = int(stringNum)
    field_vals.gridWidth = value

def buffer_val(field_vals, bufferEntry):
    stringNum = bufferEntry.get()
    value = int(stringNum)
    value = value * -1
    field_vals.bufferGiven = value

#changes the map by the user entering the entry number of
#a cell and a desired rate
def change_map(cellList, entryEntry, nitrogenEntry):
    stringNum = entryEntry.get()
    entryNum = int(stringNum)
    stringNum = nitrogenEntry.get()
    newNitrogenVal = int(stringNum)
    for cells in cellList:
        if cells.entryNum == entryNum:
            cells.nitrogen = newNitrogenVal
            break

#doesn't let the GUI run if all values haven't been entered
def complete(field_vals, file_names, root, widthEntry, heightEntry, bufferEntry):
    width_val(field_vals, widthEntry)
    height_val(field_vals, heightEntry)
    buffer_val(field_vals, bufferEntry)
    if file_names.fieldFileName != "" and file_names.yieldFileName != "" and file_names.proteinFileName != "" and file_names.choice != "" and field_vals.numYieldBins > -1 and field_vals.numProteinBins > -1 and len(field_vals.nitrogenList) > 0:
        file_names.complete = True
        root.destroy()

def closeWindow(fixWindow, file_names, satisfied):
    file_names.satisfied = satisfied
    fixWindow.destroy()

def yieldDataFrame(nameOfFile):
    #returns a dataframe of yield data with only x, y, and yield data
    yieldpoints = pandas.read_csv(nameOfFile)
    dropHeaders = [header for header in yieldpoints.columns.values if header != "x"
    and  header != "y" and  header != "yield"]
    if len(dropHeaders) != 0:
        yieldpoints = yieldpoints.drop(dropHeaders, axis=1)
    return yieldpoints

def proteinDataFrame(nameOfFile):
    #returns a dataframe of protein data with only x, y, and protein vals
    proteinpoints = pandas.read_csv(nameOfFile)
    dropHeaders = [header for header in proteinpoints.columns.values if header
    != "x" and  header != "y" and  header != "protein"]
    if len(dropHeaders) != 0:
        proteinpoints = proteinpoints.drop(dropHeaders, axis=1)
    return proteinpoints

def fieldShape(nameOfFile, bufferGiven):
    #takes all the chunks of a field and then composes them to make a giant
    #polygon out of the field boundaries
    polygons =[shape(feature['geometry']) for feature in fiona.open(nameOfFile)]
    fieldShape = cascaded_union(polygons)
    fieldShape = fieldShape.buffer(bufferGiven, cap_style = 3)
    return fieldShape

def binEqualWidth(fieldVals, yieldpoints, proteinpoints):
    #in the case of using equal bounds, this splits data in equal width parts
    #for protein and yield, saves the resulting boundaries and returns an
    #object with all those pieces of information
    fieldVals.yieldBounds = ps.esda.mapclassify.Equal_Interval(yieldpoints['yield'], k = fieldVals.numYieldBins).bins.tolist()
    fieldVals.proteinBounds = ps.esda.mapclassify.Equal_Interval(proteinpoints['protein'], k = fieldVals.numProteinBins).bins.tolist()
    fieldVals.yieldBounds[-1] = fieldVals.yieldBounds[-1] + .1
    fieldVals.proteinBounds[-1] = fieldVals.proteinBounds[-1] + .1
    return fieldVals

def binEqualDistribution(fieldVals, yieldpoints, proteinpoints):
    #splits the data into quantiles based on how many bins there are
    fieldVals.yieldBounds = ps.esda.mapclassify.Quantiles(yieldpoints['yield'], k = fieldVals.numYieldBins).bins.tolist()
    fieldVals.proteinBounds = ps.esda.mapclassify.Quantiles(proteinpoints['protein'], k = fieldVals.numProteinBins).bins.tolist()
    fieldVals.yieldBounds[-1] = fieldVals.yieldBounds[-1] + .1
    fieldVals.proteinBounds[-1] = fieldVals.proteinBounds[-1] + .1
    return fieldVals

def gridCreation(nameOfFile, gridWidth, gridHeight, field_vals):
    #This code was taken from Amy, minus me adding a few lines to it. I'm not sure
    #how this ends up working with copy right due to the fact some of it is posted online
    #it may need to be rewritten. Another possibility is doing tiles at angles which
    #is not currently implemented.

    #This function returns allthe cells created by this grid process, regardless of if they are
    #outside or not
    """
    FOLLOWING CODE THANKS TO:
    https://github.com/mlaloux/My-Python-GIS_StackExchange-answers/blob/master/Generate%20grid%20programmatically%20using%20QGIS%20from%20Python.md
    """
    fieldFilename = nameOfFile
    gridCellList = []
    xmin,ymin,xmax,ymax = fiona.open(nameOfFile).bounds #read in shapefile
    field_vals.xmin = xmin
    field_vals.ymin = ymin
    field_vals.xmax = xmax
    field_vals.ymax = ymax
    #Set grid cell dimensions
    rows = (ymax-ymin)/gridHeight
    cols = (xmax-xmin)/gridWidth
    ringXleftOrigin = xmin
    ringXrightOrigin = xmin + gridWidth
    ringYtopOrigin = ymax
    ringYbottomOrigin = ymax-gridHeight
    schema = {'geometry': 'Polygon','properties':{'fieldname':'str', 'index':'int'}}
    with fiona.open(fieldFilename+ '_grid.shp', 'w', 'ESRI Shapefile', schema) as c:
        id = 0
        for i in np.arange(cols):
            ringYtop = ringYtopOrigin
            ringYbottom =ringYbottomOrigin
            for j in np.arange(rows):
                #I added this
                newGridCell = GridCell(ringXleftOrigin, ringYbottom, ringXrightOrigin, ringYtop)
            #    newGridCell.entryNum = id
                gridCellList.append(newGridCell)
                #end of I added this
                polygon = Polygon([(ringXleftOrigin, ringYtop), (ringXrightOrigin, ringYtop), (ringXrightOrigin, ringYbottom), (ringXleftOrigin, ringYbottom)])
                c.write({'geometry':mapping(polygon),'properties':{'fieldname':fieldFilename, 'index':id}})
                ringYtop = ringYtop - gridHeight
                ringYbottom = ringYbottom - gridHeight
                id+=1
            ringXleftOrigin = ringXleftOrigin + gridWidth
            ringXrightOrigin = ringXrightOrigin + gridWidth
            id +=1
    return gridCellList

def goodCells(fieldShape, gridCellList):
    #seeing what is contained inside field
    innerCells = []
    for cells in gridCellList:
        if fieldShape.contains(cells.rectangle):
            innerCells.append(cells)
    return innerCells

#looks at the as applied map. Assumes that they're entered in order.
#this might have to be changed depending on what the shapefiles fields are
def ordering(averageCells, fileName):
    idVal = 0
    for feature in fiona.open(fileName):
        pointToFind = Point(feature['geometry']['coordinates'])
        for cells in averageCells:
            if cells.trueBounds.contains(pointToFind):
                if cells.isSorted == 0:
                    cells.entryNum = idVal
                    cells.isSorted = 1
                    idVal = idVal + 1
                break
    return averageCells

#sorts based on the entry number value
def sortCells(averageCells):
    averageCells.sort(key=operator.attrgetter('entryNum'))

def nitrogenDistribution(innerCells, bounds):
    #stratifies combinations of different bins and then assigns nitrogen rates
    newCells = []
    i = 0
    j = 0
    while i <= bounds.numYieldBins:
        while j <= bounds.numProteinBins:
            tempCellCollection = []
            for cells in innerCells:
                if cells.yieldBin == i and cells.proteinBin == j:
                    tempCellCollection.append(cells)
            newCells.extend(randomNitrogen(tempCellCollection, bounds))
            j = j + 1
        i = i + 1
        j = 0

#depending on the size of the cells and the available data, there are cases where
#protein and yield are never changed. In this case, it takes all the unbinned ones
#and goes through assigning nitrogen rates
    tempCellCollection = []
    for cells in innerCells:
        if cells.yieldBin == -1 or cells.proteinBin == -1:
            tempCellCollection.append(cells)
    newCells.extend(randomNitrogen(tempCellCollection, bounds))
    return newCells

def randomNitrogen(innerCells, bounds):
    #takes a set of cells and assigns equal-ish distribution of nitrogen
    i = 0
    random.shuffle(bounds.nitrogenList)
    for cells in innerCells:
        cells.nitrogen = bounds.nitrogenList[i % len(bounds.nitrogenList)]
        i = i + 1
    return innerCells

def averageValuesPerCell(fieldValues, innerCells, yieldpoints, proteinpoints):
    for cells in innerCells:
        #averaging yield in squares
        yieldInSquare = yieldpoints[(yieldpoints['x'] >= cells.bl_x) &
        (yieldpoints['x'] <= cells.ur_x) & (yieldpoints['y'] <= cells.ur_y) &
        (yieldpoints['y'] >= cells.bl_y)]
        mean = yieldInSquare['yield'].mean()
        cells.cropYield = mean

        #averaging protein in squares
        proteinInSquare = proteinpoints[(proteinpoints['x'] >= cells.bl_x) &
        (proteinpoints['x'] <= cells.ur_x) & (proteinpoints['y'] <= cells.ur_y) &
         (proteinpoints['y'] >= cells.bl_y)]
        mean = proteinInSquare['protein'].mean()
        cells.protein = mean

        #using predetermined yield bounds to place in yield bins
        idVal = 1
        i = 0
        while i < (len(fieldValues.yieldBounds)):
            if(cells.cropYield < fieldValues.yieldBounds[i]):
                cells.yieldBin = idVal
                break
            idVal = idVal + 1
            i = i + 1

        #using predetermined protein bounds to place in protein bins
        idVal = 1
        i = 0
        while i < (len(fieldValues.proteinBounds)):
            if(cells.protein < fieldValues.proteinBounds[i]):
                cells.proteinBin = idVal
                break
            idVal = idVal + 1
            i = i + 1

    return innerCells

def createWKTFile(fieldShape, averageCells):
    #creates a wkt with various information
    with open('OUT.csv', 'w') as csvfile:
        fieldnames = ['WKT', 'id', 'yield', 'yieldBin', 'protein', 'proteinBin', 'nitrogen']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        #writing horizontal lines
        idVal = 1
        writer.writerow({'WKT': fieldShape, 'id': idVal, 'yield': "", 'yieldBin': "",
         'protein': "", 'proteinBin': "", 'nitrogen': ""})

        # writer.writerow({'WKT': 'Baked', 'id': 'Beans', 'yield': 10})
        idVal = idVal + 1
        for cells in averageCells:
            writer.writerow({'WKT': cells.trueBounds, 'id': idVal, 'yield': cells.cropYield, 'yieldBin':
            cells.yieldBin, 'protein': cells.protein, 'proteinBin': cells.proteinBin,
             'nitrogen': cells.nitrogen})
            idVal += 1
    return

def foliumMapGenerator(file_names, field_values, innerCells):
    #generates a map based on coords
    fullField = fieldShape(file_names.fieldFileName, field_values.bufferGiven)
    xmin,ymin,xmax,ymax = fiona.open(file_names.fieldFileName).bounds
    x, y  = utm.to_latlon((xmax + xmin) / 2, (ymax + ymin) / 2, 12, 'T')
    field_map = folium.Map([x, y], zoom_start = 15)

    #applies field shape and field grid to the map
    fieldGrid = gpd.read_file(file_names.fieldFileName + "_grid.shp")
    field = gpd.read_file(file_names.fieldFileName)
    newshape = gpd.GeoDataFrame()
    newshape['geometry'] = None
    index = 0
    newshape.loc[index, 'geometry'] = transform(project, fullField)
    field_map.choropleth(geo_data = newshape.to_json(), data = None, columns = ['geometry'], line_opacity = 1, fill_opacity = 0, fill_color = "Purple")
    newdata = gpd.GeoDataFrame()
    newdata['geometry'] = None

    index = 0
    for items in fieldGrid["geometry"]:
        items = transform(project, items)
        newdata.loc[index, 'geometry'] = items
        index = index + 1

    #takes cells and makes them into something that can be applied to folium
    for cells in innerCells:
        bounds = [cells.ur_x, cells.ur_y, cells.bl_x, cells.bl_y]
        shapeOfFeature = Polygon([(cells.bl_x, cells.bl_y), (cells.ur_x, cells.bl_y), (cells.ur_x, cells.ur_y), (cells.bl_x, cells.ur_y)])
        shapeOfFeature = transform(project, shapeOfFeature)
        bl_x, bl_y, ur_x, ur_y = shapeOfFeature.bounds
        bounds = [ur_y, ur_x, bl_y, bl_x]
        field_values.nitrogenList.sort()
        fill_color = str(color(field_values, cells.nitrogen))
        folium.features.RectangleMarker(bounds = bounds, popup = "Protein: %.2f <br> Yield: %.2f <br> Nitrogen: %d <br> Entry: %d" % (cells.protein, cells.cropYield, cells.nitrogen, cells.entryNum), color = "Black", fill_color = fill_color, fill_opacity = 1).add_to(field_map)

    #popups on field. Probably want to change these
    html="""
    <p><strong> map of field </strong></p>
    """
    iframe = folium.IFrame(html=html, width=200, height=50)
    popup = folium.Popup(iframe, max_width=2650)
    x, y  = utm.to_latlon((xmax + xmin) / 2, (ymax + 25), 12, 'T')
    folium.Marker([x,y], popup=popup).add_to(field_map)
    html = """<p> Color Guide </p> <br>"""
    i = 0
    while i < len(field_values.nitrogenList):
        opening = "<p><font color = " + color(field_values, field_values.nitrogenList[i]) + ">"
        newString = "&#9608 nitrogen level of %d" % (field_values.nitrogenList[i])
        closing = "</font></p>"
        breakLine = "<br>"
        fullstring = opening + newString + closing + breakLine
        html = html + fullstring
        i = i + 1
    iframe = folium.IFrame(html=html, width=100, height=300)
    popup = folium.Popup(iframe, max_width=2650)
    x, y  = utm.to_latlon((xmin - 25) , (ymax + ymin) / 2, 12, 'T')
    folium.Marker([x,y], popup=popup).add_to(field_map)
    folium.LatLngPopup().add_to(field_map)

    #saves and opens
    field_map.save("created_map.html")
    webbrowser.open_new_tab("created_map.html")
    return

#right now, the program only supports eight max nitrogen
#values
def color(field_values, value):
    if value == field_values.nitrogenList[0]:
        return "#a10000"
    if value == field_values.nitrogenList[1]:
        return "#a15000"
    if value == field_values.nitrogenList[2]:
        return "#a1a100"
    if value == field_values.nitrogenList[3]:
        return "#416600"
    if value == field_values.nitrogenList[4]:
        return "#008282"
    if value == field_values.nitrogenList[5]:
        return "#005682"
    if value == field_values.nitrogenList[6]:
        return "#000056"
    if value == field_values.nitrogenList[7]:
        return "#6a006a"

#Just the GUI. Sets up all the classes that are
#needed to be able to run the program. CUSTOM IS
#NOT SET UP AT ALL
def guiScreen(file_names, field_vals):
    #dividing everything into frames to more easily and
    #clearly display the buttons
    newWindow = Tk()
    newWindow.title("Nitrogen Prescription Map")

    fileFrame = Frame(newWindow)
    selectionFrame = Frame(newWindow)
    logoFrame = Frame(newWindow)

    fileFrame.pack(side = TOP)
    selectionFrame.pack(side = TOP)
    logoFrame.pack(side = TOP)

    fieldFrame = Frame(fileFrame)
    yieldFrame = Frame(fileFrame)
    proteinFrame = Frame(fileFrame)
    appliedFrame = Frame(fileFrame)

    fieldFrame.pack(side = TOP)
    yieldFrame.pack(side = TOP)
    proteinFrame.pack(side = TOP)
    appliedFrame.pack(side = TOP)

    binsAndSplit = Frame(selectionFrame)
    nitrogenDimension = Frame(selectionFrame)

    binsAndSplit.pack(side = TOP)
    nitrogenDimension.pack(side = TOP)

    dataSplitting = Frame(binsAndSplit)
    yieldBins = Frame(binsAndSplit)
    proteinBins = Frame(binsAndSplit)

    dataSplitting.pack(side = LEFT)
    yieldBins.pack(side = LEFT)
    proteinBins.pack(side = LEFT)

    nitrogenInput = Frame(nitrogenDimension)
    heightWidth = Frame(nitrogenDimension)

    nitrogenInput.pack(side = LEFT)
    heightWidth.pack(side = LEFT)

    heightInput = Frame(heightWidth)
    widthInput = Frame(heightWidth)
    bufferInput = Frame(heightWidth)

    heightInput.pack(side = TOP)
    widthInput.pack(side = TOP)
    bufferInput.pack(side = TOP)

    logos = Frame(logoFrame)
    runButton = Frame(logoFrame)

    logos.pack(side = LEFT)
    runButton.pack(side = RIGHT)

    agricultureLogo = Frame(logos)
    msuLogo = Frame(logos)

    agricultureLogo.pack(side = LEFT)
    msuLogo.pack(side = LEFT)

    fieldLabel = Label(fieldFrame, text="field shape")
    fieldEntry = Entry(fieldFrame)
    fieldButton = Button(fieldFrame, text="browse", command = lambda: field_file(file_names, filedialog.askopenfilename(initialdir = "/home/neddre/Desktop/editPA", filetypes=(("???", ".shp"), ("All files","*.*"))), fieldEntry))
    fieldLabel.pack(side = LEFT)
    fieldEntry.pack(side = LEFT)
    fieldButton.pack(side = LEFT)

    yieldLabel = Label(yieldFrame, text = "yield file")
    yieldEntry = Entry(yieldFrame)
    yieldButton = Button(yieldFrame, text = "browse", command = lambda: yield_file(file_names, filedialog.askopenfilename(initialdir = "/home/neddre/Desktop/editPA", filetypes=(("???", ".csv"), ("All files","*.*"))), yieldEntry))
    yieldLabel.pack(side = LEFT)
    yieldEntry.pack(side = LEFT)
    yieldButton.pack(side = LEFT)

    proteinLabel = Label(proteinFrame, text = "protein file")
    proteinEntry = Entry(proteinFrame)
    proteinButton = Button(proteinFrame, text = "browse", command = lambda: protein_file(file_names, filedialog.askopenfilename(initialdir = "/home/neddre/Desktop/editPA", filetypes=(("???", ".csv"), ("All files","*.*"))), proteinEntry))
    proteinLabel.pack(side = LEFT)
    proteinEntry.pack(side = LEFT)
    proteinButton.pack(side = LEFT)

    appliedLabel = Label(appliedFrame, text = "applied file")
    appliedEntry = Entry(appliedFrame)
    appliedButton = Button(appliedFrame, text = "browse", command = lambda: applied_file(file_names, filedialog.askopenfilename(initialdir = "/home/neddre/Desktop/editPA", filetypes=(("???", ".shp"), ("All files","*.*"))), appliedEntry))
    appliedLabel.pack(side = LEFT)
    appliedEntry.pack(side = LEFT)
    appliedButton.pack(side = LEFT)

    options = ["equal width", "equal distribution", "custom"]
    optionVar = StringVar(dataSplitting)
    optionVar.set(options[0])
    optionBar = OptionMenu(dataSplitting, optionVar, *options, command = lambda x: option(file_names, optionVar))
    optionBar.pack(side = LEFT)

    binOptions = [1, 2, 3, 4]

    yieldBinLabel = Label(yieldBins, text = "yield bins")
    yieldBinOptionVar = IntVar(yieldBins)
    yieldBinOptionVar.set(binOptions[0])
    yieldOptions = OptionMenu(yieldBins, yieldBinOptionVar, *binOptions, command = lambda x: yieldBinSelect(field_vals, yieldBinOptionVar))
    yieldBinLabel.pack(side = TOP)
    yieldOptions.pack(side = TOP)

    proteinBinLabel = Label(proteinBins, text = "protein bins")
    proteinBinOptionVar = IntVar(proteinBins)
    proteinBinOptionVar.set(binOptions[0])
    proteinOptions = OptionMenu(proteinBins, proteinBinOptionVar, *binOptions, command = lambda x: proteinBinSelect(field_vals, proteinBinOptionVar))
    proteinBinLabel.pack(side = TOP)
    proteinOptions.pack(side = TOP)

    nitrogenLabel = Label(nitrogenInput, text = "nitrogen values")
    nitrogenEntry =  Entry(nitrogenInput)
    nitrogenSubmit = Button(nitrogenInput, text="browse", command = lambda: nitrogen_list(field_vals, nitrogenEntry, filedialog.askopenfilename(initialdir = "/home/neddre/Desktop/editPA", filetypes=(("???", ".txt"), ("All files","*.*")))))
    nitrogenLabel.pack()
    nitrogenEntry.pack()
    nitrogenSubmit.pack()

    heightLabel = Label(heightInput, text = "height (m): ")
    heightEntry = Entry(heightInput)
    heightLabel.pack()
    heightEntry.pack()

    widthLabel = Label(widthInput, text = "width (m): ")
    widthEntry = Entry(widthInput)
    widthLabel.pack()
    widthEntry.pack()

    bufferLabel = Label(bufferInput, text = "Buffer (m): ")
    bufferEntry = Entry(bufferInput)
    bufferLabel.pack()
    bufferEntry.pack()

    runProgramButton = Button(runButton, text = "run", command = lambda: complete(field_vals, file_names, newWindow, widthEntry, heightEntry, bufferEntry))
    runProgramButton.pack()

    newWindow.mainloop()

#window that lets farmers adjust cells to liking
def adjustCells(file_names, cellList):
    fixWindow = Tk()
    fixWindow.title("Fix Cells")

    changeFrame = Frame(fixWindow)
    satisfiedFrame = Frame(fixWindow)

    changeFrame.pack(side = TOP)
    satisfiedFrame.pack(side = TOP)

    entryNumFrame = Frame(changeFrame)
    newNitrogenFrame = Frame(changeFrame)
    runFrame = Frame(changeFrame)
    entryNumFrame.pack(side = TOP)
    newNitrogenFrame.pack(side = TOP)
    runFrame.pack(side = TOP)

    entryLabel = Label(entryNumFrame, text = "Entry cell: ")
    entryEntry = Entry(entryNumFrame)
    entryLabel.pack(side = LEFT)
    entryEntry.pack(side = LEFT)

    nitrogenLabel = Label(newNitrogenFrame, text = "Nitrogen value: ")
    nitrogenEntry = Entry(newNitrogenFrame)
    nitrogenLabel.pack(side = LEFT)
    nitrogenEntry.pack(side = LEFT)

    runSubmit = Button(runFrame, text = "submit", command = lambda: change_map(cellList, entryEntry, nitrogenEntry))
    runSubmit.pack()

    genMap = Frame(satisfiedFrame)
    isNowSatisfied = Frame(satisfiedFrame)

    genMap.pack(side = TOP)
    isNowSatisfied.pack(side = TOP)

    mapButton = Button(genMap, text = "generate new map", command = lambda : closeWindow(fixWindow, file_names, False))
    satisfiedButton = Button(isNowSatisfied, text = "no more changes needed", command = lambda: closeWindow(fixWindow, file_names, True))
    mapButton.pack()
    satisfiedButton.pack()

    fixWindow.mainloop()


fileNames = file_names("hello")
fieldVals = field_vals("hello")
guiScreen(fileNames, fieldVals)

#setting up the yield and protein frames
yieldframe = yieldDataFrame(fileNames.yieldFileName)
proteinframe = proteinDataFrame(fileNames.proteinFileName)
shapeOfField = fieldShape(fileNames.fieldFileName, fieldVals.bufferGiven)

#setting up bounds for yield and protein, custom is already set up
#if implemented correctly
if fileNames.choice == "equal width":
    fieldVals = binEqualWidth(fieldVals, yieldframe, proteinframe)
elif fileNames.choice == "equal distribution":
    fieldVals = binEqualDistribution(fieldVals, yieldframe, proteinframe)

#all cells are taken from the grid, next the valid ones are weeded out, means of
#yield and protein are taken
allCells = gridCreation(fileNames.fieldFileName, fieldVals.gridWidth, fieldVals.gridHeight, fieldVals)
validCells = goodCells(shapeOfField, allCells)
averageCells = averageValuesPerCell(fieldVals, validCells, yieldframe, proteinframe)

#puts all the cells in order that farmer will be applying nitrogen to them
averageCells = ordering(averageCells, fileNames.appliedFileName)
averageCells = nitrogenDistribution(averageCells, fieldVals)

createWKTFile(shapeOfField, averageCells)
foliumMapGenerator(fileNames, fieldVals, averageCells)

#lets farmer change tiles while they still want to change things
while not fileNames.satisfied:
    adjustCells(fileNames, averageCells)
    foliumMapGenerator(fileNames, fieldVals, averageCells)

exit()
