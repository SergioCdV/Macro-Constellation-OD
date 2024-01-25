%% RVHANDLE class
% This class implements the TEME Cartesian elements of a given propagated
% spacecraft
% Author: Sergio Cuevas
% Created: 2023-12-11

classdef rvhandle < handle

  properties
    r = [0 0 0].';      % Spacecraft position vector  
    v = [0 0 0].';      % Spacecraft velocity vectory
  end

end
