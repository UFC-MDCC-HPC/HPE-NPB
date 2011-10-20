using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.topology.Mesh2D;

namespace impl.common.topology.Mesh2DImpl { 

public class ICell2DImpl : BaseICell2DImpl, ICell2D
{

    public int West   { get { return X.predecessor;   } set { X.predecessor = value;   } }  // x direction
    public int East   { get { return X.successor;   } set { X.successor = value;   } }  // x direction
    
    public int North  { get { return Y.predecessor;  } set { Y.predecessor = value;  } }  // y direction
    public int South  { get { return Y.successor;  } set { Y.successor = value;  } }  // y direction
    
}

}
