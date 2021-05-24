import numpy            as np
import numpy.ctypeslib  as Flib
import ctypes, sys
import os.path

# ================================================================ #
# ===  convert__trajectoryCoordinate.py                        === #
# ================================================================ #
def convert__trajectoryCoordinate( trajectory=None, coord_xy=None ):
    # ------------------------------------------------- #
    # --- [1]   引数チェック                        --- #
    # ------------------------------------------------- #
    if ( trajectory is None ): sys.exit( "[convert__trajectoryCoordinate.py] trajectory ???" )
    if ( coord_xy   is None ): sys.exit( "[convert__trajectoryCoordinate.py] coord_xy   ???" )

    # ------------------------------------------------- #
    # --- [2]   引数準備                            --- #
    # ------------------------------------------------- #
    #  -- [2-1] 使用する引数を準備                  --  #
    ntrajectory  = int( trajectory.shape[0] )
    ncoord       = int( coord_xy.shape[0]   )
    coord_sn     = np.zeros( coord_xy.shape )
    
    #  -- [2-2] Fortranサイズへ変換                 --  #
    trajectory_  =     np.array( trajectory, dtype=np.float64    )
    coord_xy_    =     np.array( coord_xy  , dtype=np.float64    )
    coord_sn_    =     np.array( coord_sn  , dtype=np.float64    )
    ntrajectory_ = ctypes.byref( ctypes.c_int64( ntrajectory  )  )
    ncoord_      = ctypes.byref( ctypes.c_int64( ncoord       )  )

    # ------------------------------------------------- #
    # --- [3]   ライブラリをロード                  --- #
    # ------------------------------------------------- #
    #  -- [3-1] ライブラリを定義                    --  #
    pyLIB  = Flib.load_library( 'pylib.so', os.path.abspath( os.path.dirname(__file__) ) )
    
    #  -- [3-2] 入出力管理                          --  #
    pyLIB.convert__trajectorycoordinate_.argtypes = [
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        Flib.ndpointer( dtype=np.float64 ),
        ctypes.POINTER( ctypes.c_int64   ),
        ctypes.POINTER( ctypes.c_int64   ),
    ]
    pyLIB.convert__trajectorycoordinate_.restype = ctypes.c_void_p

    # ------------------------------------------------- #
    # --- [4]   関数呼出 / 返却                     --- #
    # ------------------------------------------------- #
    pyLIB.convert__trajectorycoordinate_( trajectory_, coord_xy_, coord_sn_, ntrajectory_, ncoord_ )
    return( coord_sn_ )


# ================================================================ #
# ===  テスト用 呼び出し                                       === #
# ================================================================ #
if ( __name__=='__main__' ):

    import nkUtilities.load__pointFile as lpf
    inpFile1   = "dat/coord_xy.dat"
    inpFile2   = "dat/trajectory_input.dat"
    coord_xy   = lpf.load__pointFile( inpFile=inpFile1, returnType="point" )
    trajectory = lpf.load__pointFile( inpFile=inpFile2, returnType="point" )
    
    ret        = convert__trajectoryCoordinate( trajectory=trajectory, coord_xy=coord_xy )

    outFile   = "dat/coord_sn_pyt.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=ret )

    
