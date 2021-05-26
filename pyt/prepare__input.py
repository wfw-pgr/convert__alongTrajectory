import numpy as np

# ========================================================= #
# ===  prepare__input                                   === #
# ========================================================= #

def prepare__input():

    # ------------------------------------------------- #
    # --- [1] coordinate into this                  --- #
    # ------------------------------------------------- #
    
    inpFile = "dat/extraction_ebina.dat"
    import nkUtilities.load__pointFile as lpf
    Data = lpf.load__pointFile( inpFile=inpFile, returnType="point" )

    oData      = np.zeros( ( Data.shape[0],3 ) )
    oData[:,0] = Data[:,1]
    oData[:,1] = Data[:,2]
    oData[:,2] = 0.0
    
    outFile   = "dat/trajectory_ebina.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=oData )

    
    # ------------------------------------------------- #
    # --- [2] coordinate from this                  --- #
    # ------------------------------------------------- #
    
    inpFile = "dat/extraction_terada.dat"
    import nkUtilities.load__pointFile as lpf
    Data = lpf.load__pointFile( inpFile=inpFile, returnType="point" )

    oData      = np.zeros( ( Data.shape[0],3 ) )
    oData[:,0] = Data[:,0]
    oData[:,1] = Data[:,1]
    oData[:,2] = 0.0
    
    outFile   = "dat/trajectory_terada.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=Data )


    
# ========================================================= #
# ===                                    === #
# ========================================================= #









# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    prepare__input()
