%--------------------------------------------------------------------------
% -Welcomes the user and determines whether the problem is being
%  restarted or a data file is to be read. 
% -Reads all necessary input data.
% -Initialises kinematic variables and internal variables.
%-------------------------------------------------------------------------- 
function [PRO,CON,FEM,GEOM,QUADRATURE,BC,MAT,LOAD,CONSTANT,GLOBAL,...
          PLAST,KINEMATICS,INITIAL_KINEMATICS,STRESS,Wint,Wext,Wvdamp,EmbedElt,VolumeCorrect]...
          = input_data_and_initialisation(basedir_fem,ansmlv,inputfile,CON,vtuOn)
    %--------------------------------------------------------------------------
    % Welcomes the user and determines whether the problem is being
    % restarted or a data file is to be read.
    %-------------------------------------------------------------------------- 
    PRO = welcome(basedir_fem,ansmlv,inputfile);
    fid = PRO.fid_input;
    %----------------------------------------------------------------------
    % Read input file.
    %----------------------------------------------------------------------
    [FEM,GEOM,QUADRATURE,BC,MAT,LOAD,PRO,CON,GLOBAL,EmbedElt,...
        VolumeCorrect] = reading_input_file(PRO,CON,fid);
    %----------------------------------------------------------------------
    % Obtain entities which will be constant and only computed once.
    %----------------------------------------------------------------------
    CONSTANT = constant_entities(GEOM.ndime);
    %----------------------------------------------------------------------
    % Initialise load and increment parameters.
    %----------------------------------------------------------------------
    CON.xlamb = 0;
    CON.incrm = 0;
    %----------------------------------------------------------------------
    % Initialises kinematic variables and equivalent force vector
    %---------------------------------------------------------------------
    [GEOM,LOAD,GLOBAL,PLAST,KINEMATICS,INITIAL_KINEMATICS,STRESS] = ...
         initialisation(FEM,GEOM,QUADRATURE,MAT,LOAD,CON,GLOBAL,BC,EmbedElt,VolumeCorrect);
    %----------------------------------------------------------------------
    GLOBAL.external_load_effective = GLOBAL.external_load;
    %----------------------------------------------------------------------
    %Initilize Energy totals
    %----------------------------------------------------------------------
    Wint = 0;
    Wext = 0;
    Wvdamp = 0;
    %----------------------------------------------------------------------
    % Save into restart file.
    %----------------------------------------------------------------------
    cd(PRO.job_folder);
    save_restart_file(PRO,FEM,GEOM,QUADRATURE,BC,MAT,LOAD,CON,CONSTANT,...
                      GLOBAL,PLAST,KINEMATICS,'internal')    
    if vtuOn
        output_vtu(PRO,CON,GEOM,FEM,BC,GLOBAL,MAT,PLAST,QUADRATURE,CONSTANT,KINEMATICS,INITIAL_KINEMATICS,STRESS);
    end

