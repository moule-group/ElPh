import time

def print_start():
    """ printing the start of the DCS-Flow"""
    print(" ----------------------------------------------------------------- ")
    print(" Starting on Electron Phonon Coupling (ElPh) for Organic Semiconductors !!! ")
    print(" ElPh code made by Moule Group at UC Davis (version: 0.1.2) ")
    print(" #########   $           $ ########   $        $  ")
    print(" $           $           $        $   $        $  ")
    print(" $########   $           $ #######    $########$  ")
    print(" $           $           $            $        $  ")
    print(" #########   ##########  $            $        $  ")
    print(" ----------------------------------------------------------------- ")

    print(" -------------------------- ")
    print("       *      *     *        ")
    print("                             ")
    print("   #      2#      3#         ")
    print("                             ")
    print("       *      1*    *        ")
    print("                             ")
    print("   #      #       #          ")
    print(" -------------------------- ") 
    print(time.ctime())

def print_end():
    """ printing the end of the DCS-Flow """
    print(" ----------------------------------------------------------------- ")
    print(
    """                 
    #########  ##    #  #########
    $          # #   #  $       $
    ########   #  #  #  $      $
    $          #   # #  $     $ 
    #########  #    ##  ######
    """
    )
    print(time.ctime())

def print_error(message):
    """ printing the error message 
    Args:
    message (str): The error message
    """
    print(" ----------------------------------------------------------------- ")
    print("""
    #           #
      #       #
        #   #
          #
        #   #
      #       #
    #           #
    """)
    print(message)
    print(time.ctime())