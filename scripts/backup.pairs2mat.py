import argparse

def parseArgs():
  parser = argparse.ArgumentParser(description='''This program will process read pairs
    and generate matrix of contacts.''')
  parser.add_argument('-i','--infile', dest='infile', metavar="INFILE",
                    action='store', type=str,required=True,
                    help='input file name')
  parser.add_argument('-b','--bin', dest='bin', metavar="BINSIZE",
                      action='store', type=str,required=True,
                      help='input file name')




def main()
