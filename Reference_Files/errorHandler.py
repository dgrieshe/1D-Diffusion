# Basic error handler function.  This function simply prints an error message
# and then terminates the program.

import sys

def fatalError(message):
    print('ERROR:', message.strip())
    sys.exit()
