package sp.ADI;

computation SP_ADI(x,y,z,problem_data,cells_info)
    [class = C: Class,
     instance_type = I:Instance[class = C],     
     method = MTH: SolvingMethod] 
        from adi.ADI(x,y,z,problem_data,cells_info)[class=C, instance_type=I, method=MTH]
begin
	computation txinvr:sp.solve.BlockDiagonalMatVecProduct(cells_info, problem_data)
	                             [axis=common.axis.XYZAxes, instance_type=I, class=C, method=MTH];
 
	parallel unit adi
	begin
	    slice txinvr from txinvr.matvecproduct;
	    action txinvr; 
	    protocol: seq {do copy_faces;                             
                       do compute_rhs; 
                       do txinvr; 
                       do x_solve; 
                       do y_solve; 
                       do z_solve;
                       do add}
	end


end