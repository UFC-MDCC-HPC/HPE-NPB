package adi;

application ADI_Solver3D
    [class = C: Class,
     instance_type = I:Instance[class = C],     
     method = MTH: SolvingMethod]
begin

	computation adi: ADI(x,y,z,problem_data,cells_info)[instance_type = I, class = C, method = MTH]
	computation initialize: data.Initialize(problem_data,cells_info)[instance_type = I, class = C] 
	computation lhsinit: data.LHSInit(problem_data,cells_info)[instance_type = I, class = C]
	computation exact_rhs: ExactRHS(problem_data,cells_info)[instance_type = I, class = C]
	computation verify: Verify(x,y,z,problem_data,cells_info)[instance_type = I, class = C]
	computation data_partition: MultiPartition(x,y,z,cells_info)[instance_type = I, class = C]
	computation timer: Timer

    parallel unit adi_solver3D
    begin
       slice adi from adi.adi
       slice initialize from initialize.initialize;
       slice lhsinit from lhsinit.lhsinit;
       slice exact_rhs from exact_rhs.exact_rhs;
       slice verify from verify.verify;
       slice partition from data_partition.multipartition;
       protocol: seq{ par {do initialize; do lhsinit; do partition; do exact_rhs};
                      do timer.start_timer;
                      do adi;
                      do timer.stop_timer;
                      do verify;
                      do timer.clock.show_timing
                    }
    end
 
end    
 

