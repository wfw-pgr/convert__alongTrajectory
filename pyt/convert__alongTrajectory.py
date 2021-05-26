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
    zco,nco, sco   = np.meshgrid( zcoord, ncoord, slen, indexing="ij" )
    traj_inP_sn    = np.concatenate( [sco[:,:,:,None],nco[:,:,:,None],zco[:,:,:,None]], axis=3 )
    
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
    bf = np.reshape( bf, (Lz,Ln,Ls,3) )
    b_along  = np.concatenate( [traj_inP_sn,traj_inP_xy,bf], axis=3 )
    
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


    # ------------------------------------------------- #
    # --- [5] interpolation on to new coordinate    --- #
    # ------------------------------------------------- #

    Flag3d = True
    if ( Flag3d ):
        gridData            = np.zeros( (b_along.shape[0],b_along.shape[1],b_along.shape[2],4) )
        gridData[:,:,:,0:3] = np.copy( b_along[:,:,:,0:3] )
        pointData           = np.zeros( (coord_sn.shape[0],4) )
        pointData[:,0:3]    = np.copy ( coord_sn[:,0:3] )
        field               = np.zeros( (coord_sn.shape[0],3) )
        
        for ik in range(3):
            gridData[:,:,:,3]  = np.copy( b_along[:,:,:,6+ik] )
            pointData_         = np.copy( pointData )
            ret                = itc.interpolate__tricubic( gridData=gridData, pointData=pointData_ )
            field[:,ik]        = np.copy( ret[:,3] )
            
    
    # ------------------------------------------------- #
    # --- [6] associate with xyz grid               --- #
    # ------------------------------------------------- #
    field     = np.reshape( field, (-1,3) )
    field_xyz = np.concatenate( [coord_xy,field], axis=1 )

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
