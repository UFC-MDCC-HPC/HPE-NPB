package adi;

import common.axis.*;

computation ADI(x,y,z,problem_data,cells_info) 
    [class = C: Class,
     instance_type = I:Instance[class = C],     
     method = MTH: SolvingMethod]
begin
	computation x_solve : Solver(x,problem_data,cells_info)[axis=XAxis, instance_type=I, class=C, method=MTH];
	computation y_solve : Solver(y,problem_data,cells_info)[axis=YAxis, instance_type=I, class=C, method=MTH];
	computation z_solve : Solver(z,problem_data,cells_info)[axis=ZAxis, instance_type=I, class=C, method=MTH];
	computation add: Add(problem_data,cells_info)[instance_type = I, class = C];
	computation compute_rhs: ComputeRHS(problem_data,cells_info)[instance_type = I, class = C];
	comunicator copy_faces: CopyFaces(x,y,z,problem_data,cells_info)[instance_type = I, class = C];

    parallel unit adi
    begin    
       slice x_solve from x_solve.solve;
       slice y_solve from y_solve.solve;
       slice z_solve from z_solve.solve;
       slice add from add.add
       slice compute_rhs from compute_rhs.compute_rhs
       slice copy_faces from copy_faces.copy_faces;
       action x_solve;
       action y_solve;
       action add;
       action compute_rhs; 
       action copy_faces;       
       protocol: seq {do copy_faces;                             
                      do compute_rhs; 
                      do x_solve; 
                      do y_solve; 
                      do z_solve; 
                      do add}
    end
end