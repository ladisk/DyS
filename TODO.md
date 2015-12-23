**DyS TODO list - tasks**


*Info*

- [ ] distance: positive value-separation, negative value-indentation
- [ ] If relative normal velocity of contact point is positive, than the bodies are approaching, which
corresponds to the compression phase, and if the relative normal velocity of contact point is
negative the bodies are separating, which corresponds to the restitution phase.

*Check* 

- [ ] revolute joint type: ideal (default), clearance, real, with friction (ls)
- [ ] attributes shat start with underscore are used to be modified before each simulation (probably using monte carlo simulation)
- [ ] tangent contact force not 0 when coefficient of friction is 0 (ls)
- [ ] contact point is only updated, it there is translation in contact is not detected???
- [ ] Update AABB data (AABB frame, nodes in AABB) as function of R, theta (at)
- [ ] File/Close treeview (modelview) still displays data from closed file
- [ ] Save/load with dill (or pickle)

*ToDo*

- [ ] two/multiple points body contact with penalty method (at the moment only one contact point is supported)
- [ ] before simulation starts evaluates C(q, t) vector and stops the simulation if there are non-zero elements in the array
- [ ] update display on delta time, not number of time steps (display opengl based on simulation time, not simulation step number)
- [ ] find and update (location of) new contact point between two bodies
- [ ] pause simulation
- [ ] increase time step after contact has happened
- [ ] write predictor corrector method (variable order of numerical method)
- [ ] Menu to edit (create, delete) in tree model view and automatic update in opengl widget
- [ ] Graphical representation of springs (translational, torsional)
- [ ] kinematic motors (drives) – df/dt = ?
- [ ] calculate mass moment of inertia (volume) from cloud of triangles – stl file
- [ ] use package [OOC](http://www.pythonocc.org/) for 3D visualization as it also has support for 3D modeling, neutral formats: step, iges
- [ ] visualize parameterized (with lines) geometry models
- [ ] dynamics 2D to 3D

*Done*

- [x] euler method
- [x] add dq_t - relative tangential contact velocity
- [x] Display LCS of a body
- [x] numerical error of each time step added to solution data file
- [x] contact.no_overlap() if never runs only else
- [x] Open and load project with file .dprj (Dynamic Project)
- [x] revolute joint
- [x] revolute joint to ground
- [x] fixed joint
- [x] fixed joint to ground
- [x] prismatic joint
- [x] prismatic joint to ground
- [x] spring (translational)
- [x] damping in spring (translational)
- [x] spring to ground (translational)
- [x] spring (torsional)
- [x] spring (torsional) to ground
- [x] Add bodies.txt file to read and load only bodies that are located in this file – if you have multiple occurrences of a body with same geometry you can only have one geometry file.
- [x] delete VBOs before load new project (database)
- [x] close project – empties all data from object MBD_system
- [x] load and read .dprj project file to load simulation settings (properties)
- [x] only one index buffer for all frame vertices buffer
- [x] tree model working (using model view in pyqt)
- [x] contact (bounding box and subdivide bounding box)
- [x] Save, Save As
- [x] contact (find contact, variable time step size),
- [x] Read/write text .dat files, not .txt
- [x] If load solution after finished is checked, the solution data should be loaded
- [x] Add spring energy to calculation of mechanical energy
- [x] add header data to solution file: solution_data_num.dat
- [x] contact force to zero after contact is finished
