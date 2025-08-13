function StructToVars(varsStruct)
  % Unpacks and assigns variables from structs.
  % 
  % Authors:  Ginger Schmidt (1,2), Néstor Uribe-Patarroyo (1,2) 
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA
  % 2. Institute for Medical Engineering and Science, Massachusetts Institute 
  % of Technology, 77 Massachusetts Avenue, Cambridge,, MA, USA
  % <uribepatarroyo.nestor@mgh.harvard.edu>
  
  varNames = fieldnames(varsStruct);
    
  for k = 1:numel(varNames)
    assignin('caller', varNames{k}, varsStruct.(varNames{k}));
  end
  
  
end

