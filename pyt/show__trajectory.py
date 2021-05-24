import numpy as np

# ========================================================= #
# ===  show__trajectory.py                              === #
# ========================================================= #

def show__trajectory():

    x_,y_,z_ = 0, 1, 2
    xt, yt   = 0.6, -0.7
    xc, yc   = 0.65079136497835099, -0.75925660973603493 

    # ------------------------------------------------- #
    # --- [1] load reference trajectory             --- #
    # ------------------------------------------------- #
    inpFile  = "dat/trajectory_reference.dat"
    import nkUtilities.load__pointFile as lpf
    traj_ref = lpf.load__pointFile( inpFile=inpFile, returnType="point" )

    # ------------------------------------------------- #
    # --- [2] load input trajectory                 --- #
    # ------------------------------------------------- #
    inpFile  = "dat/trajectory_input.dat"
    import nkUtilities.load__pointFile as lpf
    traj_inp = lpf.load__pointFile( inpFile=inpFile, returnType="point" )

    # ------------------------------------------------- #
    # --- [3] load coord_xy                         --- #
    # ------------------------------------------------- #
    inpFile      = "dat/bfield_ref_sn.dat"
    import nkUtilities.load__pointFile as lpf
    coord_sn     = lpf.load__pointFile( inpFile=inpFile, returnType="point" )
    
    # ------------------------------------------------- #
    # --- [4] load coord_xy                         --- #
    # ------------------------------------------------- #
    inpFile  = "dat/coord_xy.dat"
    import nkUtilities.load__pointFile as lpf
    coord_xy = lpf.load__pointFile( inpFile=inpFile, returnType="point" )
    
    # ------------------------------------------------- #
    # --- [2] plot trajectory                       --- #
    # ------------------------------------------------- #
    import nkUtilities.plot1D       as pl1
    import nkUtilities.load__config as lcf
    pngFile = "png/trajectory_input.png"
    config  = lcf.load__config()
    fig     = pl1.plot1D( config=config, pngFile=pngFile )
    fig.add__plot( xAxis=traj_ref[:,x_], yAxis=traj_ref[:,y_], \
                   label="trajectory-Ref", linestyle="--" )
    fig.add__plot( xAxis=traj_inp[:,x_], yAxis=traj_inp[:,y_], \
                   label="trajectory-Inp", linestyle="--" )
    fig.add__plot( xAxis=coord_xy[:,x_], yAxis=coord_xy[:,y_], marker=".", linewidth=0.0 )
    fig.add__plot( xAxis=coord_sn[:, 3], yAxis=coord_sn[:, 4], marker=".", linewidth=0.0, markersize=0.1 )
    # fig.add__plot( xAxis=[xt], yAxis=[yt], marker="+", label="point", linewidth=0.0 )
    # fig.add__plot( xAxis=[xc], yAxis=[yc], marker="*", label="point", linewidth=0.0 )
    fig.add__legend()
    fig.set__axis()
    fig.save__figure()


# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    show__trajectory()
