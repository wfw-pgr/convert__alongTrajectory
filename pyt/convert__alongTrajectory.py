import numpy as np
import sys, os, subprocess

# ========================================================= #
# ===  convert__alongTrajectory.py                      === #
# ========================================================= #

def convert__alongTrajectory():
    
    # ------------------------------------------------- #
    # --- [1] load bfield                           --- #
    # ------------------------------------------------- #

    inpFile = "dat/bfield_reference.dat"
    import nkUtilities.load__pointFile as lpf
    bfield  = lpf.load__pointFile( inpFile=inpFile, returnType="structured" )

    inpFile  = "dat/trajectory_reference.dat"
    import nkUtilities.load__pointFile as lpf
    traj_ref = lpf.load__pointFile( inpFile=inpFile, returnType="point" )
    
    # ------------------------------------------------- #
    # --- [2] coordinate to be interpolated         --- #
    # ------------------------------------------------- #
    #  -- [2-1] tangent vector  :: s                --  #
    ds           = np.roll( traj_ref, -1, axis=0 ) - traj_ref
    ds[-1,:]     = ds[-2,:]
    slen         = np.cumsum( np.sqrt( np.sum( ds**2, axis=1 ) ) )
    svec         = ds / np.repeat( ( np.sqrt( np.sum( ds**2, axis=1 ) ) )[:,None], 3, axis=1 )

    #  -- [2-2] normal vector   :: n                --  #
    nvec         = np.zeros( (svec.shape[0],3) )
    nvec[:,0]    = + np.copy( svec[:,1] )
    nvec[:,1]    = - np.copy( svec[:,0] )
    nvec[:,2]    =   0.0

    #  -- [2-3] make xy-coordinate of trajectory    --  #
    nc_MinMaxNum = [ -0.2, +0.2, 21 ]
    zc_MinMaxNum = [ -0.1, +0.1, 11 ]
    ncoord       = np.linspace( nc_MinMaxNum[0], nc_MinMaxNum[1], nc_MinMaxNum[2] )
    zcoord       = np.linspace( zc_MinMaxNum[0], zc_MinMaxNum[1], zc_MinMaxNum[2] )
    Ls           = traj_ref.shape[0]
    Ln           = ncoord.shape[0]
    Lz           = zc_MinMaxNum[2]
    traj_inP_xy  = np.zeros( (Lz,Ln,Ls,3) )
    for izc,zval in enumerate(zcoord):
        for inc,Anc in enumerate(ncoord):
            traj_inP_xy[izc,inc,:,:] = np.copy( traj_ref ) + Anc*nvec
        traj_inP_xy[izc,:,:,2] = zval

    #  -- [2-4] make sn-coordinate of trajectory    --  #
    # sco,nco,     = np.meshgrid( slen, ncoord, indexing="xy" )
    # print( sco.shape, nco.shape )
    zco,nco, sco   = np.meshgrid( zcoord, ncoord, slen, indexing="ij" )
    print( sco.shape, nco.shape,zco.shape )
    # print( sco.shape, nco.shape, zco.shape )
    sys.exit()
    # zco          = np.zeros_like(sco)
    traj_inP_sn  = np.concatenate( [sco[:,:,None],nco[:,:,None],zco[:,:,None]], axis=2 )
    
    # outFile      = "dat/traj_inP_xy.dat"
    # import nkUtilities.save__pointFile as spf
    # spf.save__pointFile( outFile=outFile, Data=traj_inP_xy )

    # outFile      = "dat/traj_inP_sn.dat"
    # import nkUtilities.save__pointFile as spf
    # spf.save__pointFile( outFile=outFile, Data=traj_inP_sn )
    
    # ------------------------------------------------- #
    # --- [3] interpolation                         --- #
    # ------------------------------------------------- #
    #  -- [3-1] cubic interpolation of field in xy  --  #
    import nkInterpolator.interpolate__tricubic as itc
    grid    = np.copy( bfield[:,:,:,0:3] )
    points  = np.reshape( traj_inP_xy , (-1,3) )
    bf      = np.zeros( (points.shape[0],3) )
    for ik in range(3):
        bt        = bfield[:,:,:,ik+3]
        gridData  = np.concatenate( [grid,bt[:,:,:,None]],axis=3 )
        pointData = np.concatenate( [points,np.zeros( (points.shape[0],1))], axis=1 )
        ret       = itc.interpolate__tricubic( gridData=gridData, pointData=pointData )
        bf[:,ik]  = np.copy( ret[:,3] )
        
    #  -- [3-2] combine with sn-coordinate          --  #
    bf = np.reshape( bf, (Ln,Ls,3) )
    b_along  = np.concatenate( [traj_inP_sn,traj_inP_xy,bf], axis=2 )
    
    # outFile   = "dat/bfield_ref_sn.dat"
    # import nkUtilities.save__pointFile as spf
    # spf.save__pointFile( outFile=outFile, Data=b_along )

    # ------------------------------------------------- #
    # --- [4] desired coordinate in xy              --- #
    # ------------------------------------------------- #

    #  -- [4-1] load reference trajectory           --  #
    import nkUtilities.load__pointFile as lpf
    inpFile     = "dat/trajectory_input.dat"
    trajectory  = lpf.load__pointFile( inpFile=inpFile, returnType="point" )
    
    #  -- [4-2'] load grid to be interpolated       --  #
    import nkUtilities.load__pointFile as lpf
    inpFile     = "dat/coord_xy.dat"
    coord_xy    = lpf.load__pointFile( inpFile=inpFile, returnType="point" )

    #  -- [4-3'] execute fortran program             --  #
    cmd = "./main"
    print( cmd )
    subprocess.call( cmd.split() )

    #  -- [4-4'] load converted coordinate           --  #
    inpFile     = "dat/coord_sn.dat"
    import nkUtilities.load__pointFile as lpf
    coord_sn    = lpf.load__pointFile( inpFile=inpFile, returnType="point" )

    # print( "Max, Min ( grid )" )
    # print( np.min( coord_sn[:,0] ), np.max( coord_sn[:,0] ) )
    # print( np.min( coord_sn[:,1] ), np.max( coord_sn[:,1] ) )
    
    # print( "Max, Min ( ref )" )
    # print( np.min( b_along[:,0] ), np.max( b_along[:,0] ) )
    # print( np.min( b_along[:,1] ), np.max( b_along[:,1] ) )
    
    #  -- [4-2] prepare grid to be interpolated     --  #
    # import nkUtilities.equiSpaceGrid as esg
    # x1MinMaxNum = [  0.6,  0.8,  21 ]
    # x2MinMaxNum = [ -0.8, -0.6,  21 ]
    # x3MinMaxNum = [  0.0,  0.0,   1 ]
    # coord_xy    = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
    #                                  x3MinMaxNum=x3MinMaxNum, returnType = "point" )

    # #  -- [4-3] prepare grid to be interpolated     --  #
    # import convert__trajectoryCoordinate as ctc
    # ret         = ctc.convert__trajectoryCoordinate( trajectory=trajectory, coord_xy=coord_xy )
    # print( ret.shape )

    # ------------------------------------------------- #
    # --- [5] interpolation on to new coordinate    --- #
    # ------------------------------------------------- #
    gridData            = np.zeros( (b_along.shape[0],b_along.shape[1],3) )
    gridData[:,:,0:2]   = np.copy( b_along[:,:,0:2] )
    pointData           = np.zeros( (coord_sn.shape[0],3) )
    pointData[:,0:2]    = np.copy ( coord_sn[:,0:2] )
    field               = np.zeros( (coord_sn.shape[0],3) )
    # gridData[:,:,:,0:3] = np.repeat( np.copy ( b_along[:,:,0:3] )[None,:,:,:], 5, axis=0 )
    # gridData            = np.zeros( (5,b_along.shape[0],b_along.shape[1],4) )
    # pointData           = np.zeros( (coord_sn.shape[0],4) )
    # gridData[:,:,:,0:3] = np.repeat( np.copy ( b_along[:,:,0:3] )[None,:,:,:], 5, axis=0 )
    # gridData[]
    # pointData           = np.zeros( (coord_sn.shape[0],4) )
    # pointData[:,0:3]    = np.copy ( coord_sn )
    for ik in range(3):

        # gridData[:,:,:,3] = np.copy( b_along[:,:,6+ik] )
        gridData[:,:,2]   = np.copy( b_along[:,:,6+ik] )
        pointData_        = np.copy( pointData )

        outFile   = "dat/gridData.dat"
        import nkUtilities.save__pointFile as spf
        spf.save__pointFile( outFile=outFile, Data=gridData )

        outFile   = "dat/pointData.dat"
        import nkUtilities.save__pointFile as spf
        spf.save__pointFile( outFile=outFile, Data=pointData_ )

        # print( np.min( gridData[:,:,:,3] ), np.max( gridData[:,:,:,3] )  )
        print( gridData.shape )
        print( pointData_.shape )
        # print( np.min( gridData[:,:,:,0] ), np.max( gridData[:,:,:,0] )  )
        # print( np.min( gridData[:,:,:,1] ), np.max( gridData[:,:,:,1] )  )
        # print( np.min( pointData_[:,0] ), np.max( pointData_[:,0] )  )
        # print( np.min( pointData_[:,1] ), np.max( pointData_[:,1] )  )
        import nkInterpolator.interpolate__bilinear as bil
        ret = bil.interpolate__bilinear( gridData=gridData, pointData=pointData_ )
        # ret               = itc.interpolate__tricubic( gridData=gridData, pointData=pointData_ )
        field[:,ik]       = np.copy( ret[:,2] )
        print( np.min( ret[:,2] ), np.max( ret[:,2] ) )
    
    # ------------------------------------------------- #
    # --- [6] associate with xyz grid               --- #
    # ------------------------------------------------- #
    field     = np.reshape( field, (-1,3) )
    field_xyz = np.concatenate( [coord_xy,field], axis=1 )
    # print( coord_xy.shape, field.shape, field_snz.shape )
    # print( np.min( field[:,0] ), np.max( field[:,0] ) )
    # print( np.min( field[:,1] ), np.max( field[:,1] ) )
    # print( np.min( field[:,2] ), np.max( field[:,2] ) )

    outFile   = "dat/field_xyz.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=field_xyz )
    
    print( field_xyz.shape )
    return()

# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    convert__alongTrajectory()
