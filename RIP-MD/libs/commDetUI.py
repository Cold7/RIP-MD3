from functools import reduce
import readline # Allows the correct use of arrow keys with the interactive interface.
from os import name as osname
from subprocess import run as sprun
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import use as mpluse
from scipy.cluster import hierarchy
import time

import commDet



# Exception to warn if a command with errors was entered.
class BadCommand(Exception):
    def __init__(self,comando,msg):
        self.comando = comando # Faulty command.
        self.msg = msg # Message to the user.


# Dictionary to translate abbreviations to the interaction labels in the graphs.
interactionDict = {"ca":"Alpha Carbon Distance","s":"Salt Bridge","rr":"Arg-Arg","H":"Hydrogen Bond","ss":"Disulfide Bridge","cp":"Cation-Pi","pp":"Pi-Pi","c":"Coulomb","v":"Van der Waals"}


##################################################################
##
## Parses optional aguments in commands. Returns argDict.
##
## argDict: When an option is entered with no arguments, it
## will have a single space. If it had arguments, it will
## start with two spaces and then it will have the arguments.
## Else it will be empty.
##
##################################################################
def parseOptions(inputList,options,comando):
    # Initialize argDict
    argDict = {}
    for option in options:
        argDict[option] = ""

    # Parse input
    opt = None
    i = 1
    while i < len(inputList):
        try:
            while inputList[i] not in argDict.keys():
                argDict[opt] += " " + inputList[i]
                i += 1
        except IndexError:
            break
        except KeyError:
            raise BadCommand(comando, "argument '" + inputList[i] + "' not recognized.")

        opt = inputList[i]
        if argDict[opt]:
            raise BadCommand(comando, "option " + opt + " has been entered more than once.")
        argDict[opt] += " "  # Added to know which options were entered.
        i += 1
    return argDict

####################################################
##
## Main function for the interactive CLI.
##
####################################################
def cdUI(outputFolder,backend):
    print("Welcome to RIP-MD community detection analysis.\nUse command help for more information.\n")

    # Setting the matplotlib backend to default.
    try:
        mpluse(backend,force=True)
    except ImportError:
        print("Warning: Matplotlib couldn't setup the backend properly, show plot capabilities unavailable.")
        try:
            mpluse("agg",force=True)
        except ImportError:
            print("Error: problem persists... exiting.\n")
            exit()
    except ValueError:
        print("Warning: Matplotlib couldn't setup the backend properly, show plot capabilities unavailable.")
        try:
            mpluse("agg", force=True)
        except ImportError:
            print("Error: problem persists... exiting.\n")
            exit()

    # Checking for consensus.gml in the output folder.
    if outputFolder:
        try:
            graphPath = outputFolder + "/RIP-MD_Results/Graphs/consensus.gml" # Now RIP-MD always saves a consensus.gml file for this purpose.
            grafo = nx.read_gml(graphPath)
            print("Loaded graph from " + graphPath)
        except FileNotFoundError:
            graphPath = None
            print("Warning: couldn't find the consensus.gml file in the output folder.")
        except:
            graphPath = None
            print("Warning: couldn't load the provided GML file.")
    else:
        print("Warning: no output folder specified, you will need to load a GML file manually.")
        graphPath = None

    print("")

    # Main loop for UI.
    CDhistory = open("CDhistory.log", "a") # To log information on saved files.
    while True:
        try:
            ui = input("RIP-MD>> ")  # User input.
            for command in ui.split(";"): # ";" allows executing multiple commands secuentally in a single line. # Example: "load consensus.gml ; girvanNewman;plot -c colorful -t 3".
                if command=="":
                    continue

                # Standardizing spaces.
                t = 1
                pui = ""
                for car in command:
                    if car == " ":
                        t = 1
                    elif t == 1:
                        pui += " " + car
                        t = 0
                    else:
                        pui += car
                pui = pui[1:].split(" ") # Parsed user input list.

                # Exits the UI.
                if pui[0] == "exit": # exit
                    print(" ".join(pui[1:]))
                    exit()

                # Clears the screen. Should work in most OSs.
                elif len(pui) == 1 and pui[0] == "clear": # clear
                    if osname.lower() == "posix":
                        sprun("clear")
                    elif osname.lower() == "nt":
                        sprun("cls")
                    else:
                        pass

                # Shows information about the module, its commands, their syntax and examples.
                elif pui[0] == "help": # help
                    pass # TODO

                # Loads a specified graph in GML format. Uses GML as default since it keeps node and edge data.
                elif pui[0] == "load": # load [path]
                    if len(pui) == 1:
                        try:
                            graphPath = outputFolder + "/RIP-MD_Results/Graphs/consensus.gml" # Now RIP-MD always saves a consensus.gml file for this purpose.
                            grafo = nx.read_gml(graphPath)
                            print("Loaded graph from " + graphPath)
                            continue # fixed after sending final
                        except TypeError:
                            raise BadCommand(command,"no output folder specified at startup.")
                        except FileNotFoundError:
                            graphPath = None
                            raise BadCommand(command,"couldn't find the consensus.gml file in the output folder.")
                        except:
                            graphPath = None
                            raise BadCommand(command,"couldn't load the provided GML file.")

                    graphPath = " ".join(pui[1:])
                    try:
                        grafo = nx.read_gml(graphPath)
                        print("Loaded graph from " + graphPath)
                    except FileNotFoundError:
                        graphPath = None
                        raise BadCommand(command,"couldn't find the provided GML file.")
                    except:
                        graphPath = None
                        raise BadCommand(command,"couldn't load the provided GML file.")

                # Prints the path of the current loaded graph.
                elif pui[0] == "where": # where
                    if graphPath:
                        print("Graph loaded from " + graphPath)
                    else:
                        print("No graph loaded.")

                # Runs girvanNewman algorithm on the loaded graph and saves etrDict and nodes. -k option finds best graph with commNumber communities. -b option finds graph with the highest modularity.
                elif pui[0] == "girvanNewman": # girvanNewman [-i ca s rr H ss cp pp c v][w] [-t tiempo] # girvanNewman -k commNumber [-s path] # girvanNewman -b [-s path] # Each interaction can be specified a positive weight.
                    argDict = parseOptions(pui, {"-i", "-t","-k","-b","-s"},command)

                    if (argDict["-i"] or argDict["-t"]) and (argDict["-k"] or argDict["-b"]):
                        raise BadCommand(command,"incorrect mix of options.\n-i and -t options are for running the algorithm.\n-b and -k are for analysing the algorithm's output.")
                    if argDict["-k"] and argDict["-b"]:
                        raise BadCommand(command,"-k and -b options have to be used separately.")
                    if argDict["-s"] and not (argDict["-k"] or argDict["-b"]):
                        raise BadCommand(command, "-s option only works with -b or -k.")

                    # Parsing interactions
                    weightInt = None
                    interactions = dict() # (interactionType,weight)
                    if argDict["-i"]:
                        try:
                            iterargs=iter(argDict["-i"][2:].split(" "))
                            intarg = next(iterargs)
                            while True:
                                if interactionDict[intarg] in interactions.keys():
                                    raise BadCommand(command,"interaction "+ intarg +" entered more than once.")
                                interactions[interactionDict[intarg]] = 1
                                nintarg = next(iterargs)
                                if nintarg not in interactionDict.keys():
                                    weight = float(nintarg)
                                    if weight <= 0:
                                        raise ValueError
                                    interactions[interactionDict[intarg]] = weight
                                    intarg = next(iterargs)
                                    weightInt=True
                                else:
                                    intarg = nintarg
                        except StopIteration:
                            pass
                        except KeyError:
                            raise BadCommand(command, "interaction type '" + intarg + "' not recognized.")
                        except ValueError:
                            raise BadCommand(command, "weights should be strictly positive numbers.")
                    weightInt = dict(interactions) if weightInt else None

                    # Setting time threshold
                    try:
                        tiempo = float(argDict["-t"][2:]) if  argDict["-t"] else 75
                        if not (0 <= tiempo <= 100):
                            raise ValueError
                    except ValueError:
                        raise BadCommand(command, "time threshold should be a float between 0 and 100.")

                    # Boolean for finding community structure with the best modularity.
                    findBestMod = bool(argDict["-b"])
                    if len(argDict["-b"]) > 1:
                        raise BadCommand(command, "option -b expects no arguments.")

                    # Setting option to find k communities
                    try:
                        k = int(argDict["-k"]) if argDict["-k"] else None
                        if k and k <= 0:
                            raise ValueError
                    except ValueError:
                        raise BadCommand(command,"option -k requires a positive integer.")

                    # Setting path for saving graph
                    bPath = argDict["-s"][2:] if argDict["-s"][2:] else "GNbestModGraph.gml"
                    bPath = bPath + ".gml" if bPath[-4:] != ".gml" else bPath

                    kPath = argDict["-s"][2:] if argDict["-s"][2:] else "GNkCommGraph.gml"
                    kPath = kPath + ".gml" if kPath[-4:] != ".gml" else kPath

                    # Running the algorithm
                    if not (k or findBestMod):
                        try:
                            print("Started running Girvan-Newman algorithm. Time: "+str(time.strftime("%H:%M:%S"))+".")
                            etrDict, nodes = commDet.girvanNewman(grafo, interactions.keys(), time=tiempo,weights=weightInt)
                            if len(etrDict) == 0:
                                del etrDict
                                del nodes
                                raise BadCommand(command,"your filters have rendered an empty graph.")
                            else:
                                print("Girvan-Newman algorithm run correctly. Time: "+str(time.strftime("%H:%M:%S"))+".")
                                currentGirvanNewman = (graphPath," ".join(pui))
                        except NameError:
                            raise BadCommand(command,"you need to load a graph before running an algorithm.")

                    # Finding community structure as required
                    if findBestMod:
                        try:
                            G, mod,commNumber = commDet.findBestModGN(etrDict, nodes)
                            print("Found division in "+ str(commNumber) +" communities with best modularity: " + str(mod))

                            nx.write_gml(G,bPath)
                            print("Community labeled graph saved to '" + bPath + "'")
                            CDhistory.write("\nDate: " + str(time.strftime("%y/%m/%d")) + " Time: " + str(time.strftime("%H:%M:%S")) + "\n" +
                                        "With graph loaded from '" + currentGirvanNewman[0] + "'\n" +
                                        "and girvanNewman run as '" + currentGirvanNewman[1] + "'\n" +
                                        "Found division in " + str(commNumber) + " communities with best modularity: " + str(mod) + "\n"
                                        "Community labeled graph saved to '" + bPath + "'\n")
                        except FileNotFoundError:
                            print("Warning: couldn't save best modularity graph. No such file or directory: '" + bPath + "'.")
                        except NameError:
                            raise BadCommand(command,"you need to first run girvanNewman algorithm and then run the -b option.")
                    if k:
                        try:
                            G,outOfBounds,final_k = commDet.findKCommGN(etrDict, nodes, k)

                            if outOfBounds=="cc":
                                print("Warning: k smaller than number of connected components of filtered graph: " + str(final_k))
                            elif outOfBounds=="nodes":
                                print("Warning: k bigger than number of nodes in filtered graph: " + str(final_k))

                            nx.write_gml(G,kPath)
                            print("Community labeled graph saved to '" + kPath + "'. Modularity: " + str(commDet.modularity(G)))
                            CDhistory.write("\nDate: " + str(time.strftime("%y/%m/%d")) + " Time: " + str(time.strftime("%H:%M:%S")) + "\n" +
                                        "With graph loaded from '" + currentGirvanNewman[0] + "'\n" +
                                        "girvanNewman run as '" + currentGirvanNewman[1] + "'\n" +
                                        "Found division in " + str(final_k) +" communities. Modularity: " + str(commDet.modularity(G)) + "\n" +
                                        "Community labeled graph saved to '" + kPath + "'\n")
                            if outOfBounds=="cc":
                                CDhistory.write("Warning: k smaller than number of connected components of filtered graph: " + str(final_k) + "\n")
                            elif outOfBounds=="nodes":
                                CDhistory.write("Warning: k bigger than number of nodes in filtered graph: " + str(final_k) + "\n")
                        except FileNotFoundError:
                            print("Warning: couldn't save k communities graph. No such file or directory: '" + kPath + "'")
                        except NameError:
                            raise BadCommand(command,"you need to first run girvanNewman algorithm and then run the -k option.")

                # Plots a dendrogram with the output of girvanNewman.
                elif pui[0] == "dendrogram": # dendrogram [-n] [-s [savepath]] [-j jumpfunc] [-c coloropt] [-t thresholdOrcommnumber] [-d dpi] [-z xsize,ysize] [-nb]
                    argDict = parseOptions(pui,{"-n","-j","-c","-s","-t","-d","-z","-nb"},command)

                    # Boolean values for showing or saving the graph
                    show=not bool(argDict["-n"])
                    if len(argDict["-n"]) > 1:
                        raise BadCommand(command,"option -n expects no arguments.")
                    save=bool(argDict["-s"])
                    if not (show or save):
                        raise BadCommand(command,"you must either show or save the dendrogram.")

                    # Setting the threshold
                    try:
                        umbral = float(argDict["-t"]) if argDict["-t"] else 0.7
                        if umbral >= 1:
                            umbral = int(umbral)
                        elif umbral < 0:
                            raise ValueError
                    except ValueError:
                        raise BadCommand(command,"threshold should be either an integer or a float between 0 and 1.")

                    # Setting the DPI
                    try:
                        DPI = float(argDict["-d"]) if argDict["-d"] else 300.0
                        if DPI <= 0:
                            raise ValueError
                    except ValueError:
                        raise BadCommand(command,"DPI should be a positive number.")

                    # Setting the color option
                    copt = argDict["-c"][2:] if argDict["-c"] else None
                    if not (copt is None or copt in {"simple", "colorful", "alternate"}):
                        raise BadCommand(command,"unrecognized option: '-c "+argDict["-c"]+"'.")
                    if copt == "none":
                        copt = None
                    # Setting outside plot color transparent if required
                    outColor = "#00000000" if argDict["-nb"] else "white"
                    if len(argDict["-nb"]) > 1:
                        raise BadCommand(command, "option -nb expects no arguments.")

                    # Setting the path for saving the dendrogram
                    dendPath = argDict["-s"][2:] if argDict["-s"][2:] else "dendrogram.svg"

                    # Setting the figure size
                    try:
                        size = tuple(map(float,argDict["-z"].split(","))) if argDict["-z"] else (16.0,9.0)
                        if not (len(size)==2 and size[0] > 0  and size[1] > 0):
                            raise ValueError
                    except ValueError:
                        raise BadCommand(command,"figure size should be two positive floats separated by a comma.")

                    # Setting the jump function
                    if argDict["-j"]:
                        if argDict["-j"][2:] == "constant":
                            jumpf=lambda x: 1
                        elif argDict["-j"][2:] == "sqrt":
                            jumpf = lambda x: np.sqrt(x)
                        elif argDict["-j"][2:] == "linear":
                            jumpf=lambda x: x
                        else:
                            raise BadCommand(command,"jump function not recognized.")
                    else:
                        jumpf = lambda x: np.sqrt(x)

                    try:
                        Z, leafLabels, colors, stopColoring = commDet.setupDendrogram(etrDict, nodes, jump=jumpf, colored=copt, threshold=umbral)
                    except NameError:
                        raise BadCommand(command,"you need to run girvanNewman before plotting a dendrogram.")

                    plt.figure(figsize=size, facecolor=outColor)
                    hierarchy.dendrogram(Z, link_color_func=lambda u: colors[u] if u < stopColoring else '#00000000', labels=leafLabels)

                    if save:
                        try:
                            plt.savefig(dendPath, dpi=DPI)
                            print("Dendrogram saved as '" + dendPath + "'")
                            CDhistory.write("\nDate: " + str(time.strftime("%y/%m/%d")) + " Time: " + str(time.strftime("%H:%M:%S")) + "\n" +
                                        "With graph loaded from '" + currentGirvanNewman[0] + "'\n" +
                                        "girvanNewman run as '" + currentGirvanNewman[1] + "'\n" +
                                        "dendrogram run as '" + command + "'\n"
                                        "Dendrogram saved as '" + dendPath + "'\n")
                        except FileNotFoundError:
                            print("Warning: couldn't save figure. No such file or directory: '" + dendPath + "'.")
                        except:
                            raise BadCommand(command,"couldn't save figure.") # fixed after sending final
                    if show:
                        plt.show()

                # Runs the labelPropagation algorithm on the loaded graph and saves the resulting graph in GML format.
                elif pui[0] == "labelProp": # labelProp [-i ca s rr H ss cp pp c v][w] [-t time] [-c setOfNodes] [-f] [-s path] [-n] [-r] # Each interaction can be specified a positive weight.
                    argDict = parseOptions(pui,{"-i", "-t", "-c", "-f","-s","-n","-r"},command)

                    # Parsing interactions
                    weightInt = None
                    interactions = dict()  # (interactionType,weight)
                    if argDict["-i"]:
                        try:
                            iterargs = iter(argDict["-i"][2:].split(" "))
                            intarg = next(iterargs)
                            while True:
                                if interactionDict[intarg] in interactions.keys():
                                    raise BadCommand(command, "interaction " + intarg + " entered more than once.")
                                interactions[interactionDict[intarg]] = 1
                                nintarg = next(iterargs)
                                if nintarg not in interactionDict.keys():
                                    weight = float(nintarg)
                                    if weight <= 0:
                                        raise ValueError
                                    interactions[interactionDict[intarg]] = weight
                                    intarg = next(iterargs)
                                    weightInt = True
                                else:
                                    intarg = nintarg
                        except StopIteration:
                            pass
                        except KeyError:
                            raise BadCommand(command, "interaction type '" + intarg + "' not recognized.")
                        except ValueError:
                            raise BadCommand(command, "weights should be strictly positive numbers.")
                    weightInt = dict(interactions) if weightInt else None

                    # Setting time threshold
                    try:
                        tiempo = float(argDict["-t"][2:]) if  argDict["-t"] else 75
                        if not (0 <= tiempo <= 100):
                            raise ValueError
                    except ValueError:
                        raise BadCommand(command, "time threshold should be a float between 0 and 100.")

                    # Boolean for using an enumeration instead of node names as labels.
                    nal = not bool(argDict["-n"])
                    if len(argDict["-n"]) > 1:
                        raise BadCommand(command,"option -n expects no arguments.")

                    # Boolean for removing isolated nodes.
                    ran = bool(argDict["-r"])
                    if len(argDict["-r"]) > 1:
                        raise BadCommand(command,"option -r expects no arguments.")

                    # Boolean for fixing labels of centNodes.
                    force = bool(argDict["-f"])
                    if len(argDict["-f"]) > 1:
                        raise BadCommand(command,"option -f expects no arguments.")
                    if force and not argDict["-c"]:
                        raise BadCommand(command, "option -f only works if -c is used.")

                    # Setting the path for saving the graph
                    labelPropPath = argDict["-s"][2:] if argDict["-s"][2:] else "labelPropGraph.gml"
                    labelPropPath = labelPropPath + ".gml" if labelPropPath[-4:] != ".gml" else labelPropPath

                    # Parsing -c option
                    centNodes = set()
                    if argDict["-c"]: # Node names shouldn't have any character in the set {" ","\"","'","{","}","(",")","[","]",","}
                        parsedc=""
                        for car in argDict["-c"][2:]:
                            if car not in {" ","\"","'","{","}","(",")","[","]"}:
                                parsedc += car
                        centNodes = set(parsedc.split(","))

                    # Running the algorithm
                    try:
                        print("Started running label propagation algorithm. Time: " + str(time.strftime("%H:%M:%S"))+".")
                        labelPropG = commDet.labelPropagation(grafo, interactions.keys(), tiempo, centNodes, force, nal, ran,weights=weightInt)
                        print("Label propagation algorithm run correctly. Time: " + str(time.strftime("%H:%M:%S"))+".")
                        nx.write_gml(labelPropG, labelPropPath)
                        print("Community labeled graph saved to '" + labelPropPath +"'. Modularity: " + str(commDet.modularity(labelPropG)))
                        CDhistory.write("\nDate: " + str(time.strftime("%y/%m/%d")) + " Time: " + str(time.strftime("%H:%M:%S")) + "\n" +
                                        "With graph loaded from '" + graphPath + "'\n" +
                                        "labelprop run as '" + " ".join(pui) + "'\n" +
                                        "Community labeled graph saved to '" + labelPropPath +"'. Modularity: " + str(commDet.modularity(labelPropG)) + "'\n")
                    except NameError:
                        raise BadCommand(command,"you need to load a graph before running an algorithm.")
                    except FileNotFoundError:
                        print("Warning: couldn't save the returned graph. No such file or directory: '" + labelPropPath + "'.")

                else:
                    raise BadCommand(command,"command not found.")
        except BadCommand as bcex:
            print(bcex.comando + ": " + bcex.msg)
            continue
        except KeyboardInterrupt:
            print("\nKeyboardInterrupt")
            continue
        except EOFError:
            CDhistory.close()
            print("\n")
            exit()
        except SystemExit:
            CDhistory.close()
            raise
        except:
            CDhistory.close()
            raise
