import numpy as np

# ========================================================= #
# ===  show__trajectory.py                              === #
# ========================================================= #

def show__trajectory():

    x_,y_,z_ = 0, 1, 2
    xt, yt   = 0.6, -0.7

    # ------------------------------------------------- #
    # --- [1] load trajectory                       --- #
    # ------------------------------------------------- #
    inpFile = "dat/trajectory_input.dat"
    import nkUtilities.load__pointFile as lpf
    Data = lpf.load__pointFile( inpFile=inpFile, returnType="point" )

    # ------------------------------------------------- #
    # --- [2] plot trajectory                       --- #
    # ------------------------------------------------- #
    import nkUtilities.plot1D       as pl1
    import nkUtilities.load__config as lcf
    pngFile = "png/trajectory_input.png"
    config  = lcf.load__config()
    fig     = pl1.plot1D( config=config, pngFile=pngFile )
    fig.add__plot( xAxis=Data[:,x_], yAxis=Data[:,y_], label="trajectory" )
    fig.add__plot( xAxis=[xt], yAxis=[yt], marker="+", label="point", linewidth=0.0 )
    fig.add__legend()
    fig.set__axis()
    fig.save__figure()


# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    show__trajectory()
