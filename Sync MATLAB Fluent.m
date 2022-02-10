fclose('all');
clear
clc
%% defining properties

g=9.81; Th=315; Tc=305; eps=0.7; phi=0.03;
rhof=993; kf=0.628; cpf=4178; betaf=36.2*10^(-5); muf=695*10^(-6); df=0.385; alphaf=kf/(rhof*cpf);
rhop=3970; kp=36; cpp=765; betap=0.85*10^(-5); dp=25; alphap=kp/(rhop*cpp);
result=zeros(3,4,5,5,2); Nusselt=zeros(3,4,5,5,2); s0=zeros(3,4,5,5,2); s1=zeros(3,4,5,5,2); stot=zeros(3,4,5,5,2); c=zeros(3,5); perm=zeros(3,4);

Ra=[10^4,10^5,10^6];             %rayleigh number
Da=[10^(-4),0.001,0.01,0.1];      %Darcy number
n=[0.6,0.8,1,1.2,1.4];          %power-law index
prsk=[0,0.25,0.5,0.75,1];        %partially filled factor
k_potop=[1.1162,0.2440];                  %ratio of kporous to kp
Pr=10;                           %Prandtl number

R=((Ra*alphaf^2)/(g*betaf*(Th-Tc))).^(1/3); %Cavity radius
for i=1:3
    perm(i,:)=R(i)^2.*Da;  %permeability, each row for one radius
end
viscous_resistance=1./perm;  %viscous resistance, each row for one radius

for i=1:3
     for j=1:5
        c(i,j)=(Pr*rhof*R(i)^(2*n(j)-2))/alphaf^(n(j)-2); %consistency index, each row for one diameter & column for n
     end
end
d=((perm*150*(1-eps)^2)/eps^3).^0.5; %each row for one radius
c2=(3.5*(1-eps))./(d*eps^3);   %Inertial coeff, each row for one radius

%% Executing cases
for a1=1:3                              %Rayleigh/Radius
    for a2=1:4                          %Darcy
        for a3=1:5                      %n
            for a4=2:4                  %partially filled factor
                for a5=1:2              %ratio of kporous to kp
                    %% modifying Fluent UDFs
%                     property('propUDF.c',df,dp,n(a3),c(a1,a3),rhof,rhop,kp,kf,cpf); %creating Fluent property udf for each case
%                     entropy('entropyUDF.c',Th,Tc,R(a1),kf); %creating Fluent entropy udf for each case
                    %% modifying Fluent journal
                    jou = fopen('j.jou','w'); %create/open journal file
                    fprintf(jou,'(cx-gui-do cx-activate-item "MenuBar*FunctionsSubMenu*Compiled...")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-text-entry "Compiled UDFs*TextEntry2(Library Name)" "libudf%s%s")',num2str(a1),num2str(a3));%reading Compiled properties UDF
%                     fprintf(jou,'\r\n(cx-gui-do cx-add-list-items "Select File*List" ''("propUDF.c") #f)'); 
%                     fprintf(jou,'\r\n(cx-gui-do cx-set-text-entry "Select File*FilterText" "\\*")');
%                     fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Select File*Apply")');
%                     fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Select File*OK")');
%                     fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Compiled UDFs*PushButton3(Build)")');
%                     fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Question*OK")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Compiled UDFs*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "MenuBar*FunctionsSubMenu*Interpreted...")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Interpreted UDFs*PushButton3(Browse)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-text-entry "Select File*Text" "entropyUDF%s.c")',num2str(a1)); %Interpreting entropy UDF
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Select File*OK")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Interpreted UDFs*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Interpreted UDFs*PanelButtons*PushButton2(Cancel)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "MenuBar*ReadSubMenu*Mesh...")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-text-entry "Select File*FilterText" "C:\\Users\\---\\mesh\\*")'); %mesh directory
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Select File*Apply")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-text-entry "Select File*Text" "%s.msh")',num2str(prsk(a4))); %read different meshes for each a4
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Select File*OK")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "General*Frame1*Table1*Frame1(Mesh)*ButtonBox1(Mesh)*PushButton1(Scale)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-toggle-button "Scale Mesh*Frame2(Scaling)*Table2(Scaling)*ToggleBox1*Specify Scaling Factors" #f)');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Scale Mesh*Frame2(Scaling)*Table2(Scaling)*ToggleBox1*Specify Scaling Factors")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Scale Mesh*Frame2(Scaling)*Table2(Scaling)*Frame3(Scaling Factors)*RealEntry1(X)" ''( %s))',num2str(R(a1)));
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Scale Mesh*Frame2(Scaling)*Table2(Scaling)*Frame3(Scaling Factors)*RealEntry2(Y)" ''( %s))',num2str(R(a1)));
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Scale Mesh*Frame2(Scaling)*Table2(Scaling)*PushButton4(Scale)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Scale Mesh*PanelButtons*PushButton1(Close)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-toggle-button "General*Frame1*Table1*Frame3*Frame1*CheckButton1(Gravity)" #f)');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "General*Frame1*Table1*Frame3*Frame1*CheckButton1(Gravity)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "General*Frame1*Table1*Frame3*Frame1*Frame2(Gravitational Acceleration)*RealEntry2(Y)" ''( -%s))',num2str(g));
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "General*Frame1*Table1*Frame3*Frame1*Frame2(Gravitational Acceleration)*RealEntry2(Y)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Setup|Models"))');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Models*Frame1*Table1*Frame1*List1(Models)" ''( 0))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Models*Frame1*Table1*Frame1*List1(Models)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Models*Frame1*Table1*PushButton2(Edit)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-toggle-button "Multiphase Model*Frame1*Table1*Frame1(Model)*ToggleBox1(Model)*Volume of Fluid" #f)');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Multiphase Model*Frame1*Table1*Frame1(Model)*ToggleBox1(Model)*Volume of Fluid")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-position "Multiphase Model" ''(x 290 y 255))');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-position "Multiphase Model" ''(x 290 y 209))');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-toggle-button "Multiphase Model*Frame1*Table1*Frame1(Model)*ToggleBox1(Model)*Mixture" #f)');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Multiphase Model*Frame1*Table1*Frame1(Model)*ToggleBox1(Model)*Mixture")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Multiphase Model*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Models*Frame1*Table1*Frame1*List1(Models)" ''( 1))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Models*Frame1*Table1*Frame1*List1(Models)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Models*Frame1*Table1*PushButton2(Edit)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-toggle-button "Energy*Frame1(Energy)*Table1(Energy)*Frame1*ToggleBox1*CheckButton1(Energy Equation)" #f)');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Energy*Frame1(Energy)*Table1(Energy)*Frame1*ToggleBox1*CheckButton1(Energy Equation)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Energy*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Information*OK")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Setup|Materials"))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Materials*Frame1*Table1*ButtonBox2*PushButton1(Create/Edit)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Create/Edit Materials*Frame1*Table1*Frame1*ButtonBox3*PushButton2(User-Defined Database)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Open Database*PanelButtons*PushButton2(Cancel)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Create/Edit Materials*Frame1*Table1*Frame1*ButtonBox3*PushButton1(Fluent Database)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Database Materials*Frame1*Table1*Frame1*Frame1*List1(Materials)" ''( 559))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Database Materials*Frame1*Table1*Frame1*Frame1*List1(Materials)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-position "Database Materials" ''(x 290 y 218))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Database Materials*PanelButtons*PushButton1(Copy)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Database Materials*PanelButtons*PushButton1(Close)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame4*Frame2*DropDownList1" ''( 7))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame4*Frame2*DropDownList1")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame4*Frame2*RealEntry3" ''( %s))',num2str(rhof));
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame5*Frame2*RealEntry3" ''( %s))',num2str(cpf));
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame6*Frame2*DropDownList1" ''( 5))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame6*Frame2*DropDownList1")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "User-Defined Functions*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame7*Frame2*DropDownList1" ''( 11))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame7*Frame2*DropDownList1")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "User-Defined Functions*Frame1*List1" ''( 1))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "User-Defined Functions*Frame1*List1")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "User-Defined Functions*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame21*Frame2*RealEntry3" ''( %s))',num2str(betaf));
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Create/Edit Materials*PanelButtons*PushButton1(Change/Create)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Create/Edit Materials*PanelButtons*PushButton1(Close)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Materials*Frame1*Table1*ButtonBox2*PushButton1(Create/Edit)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Create/Edit Materials*PanelButtons*PushButton1(Close)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Materials*Frame1*Table1*Frame1*List1(Materials)" ''( 2))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Materials*Frame1*Table1*Frame1*List1(Materials)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Materials*Frame1*Table1*ButtonBox2*PushButton1(Create/Edit)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-text-entry "Create/Edit Materials*Frame1*Table1*Frame1*Frame1*Table1*TextEntry1(Name)" "al2o3")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame4*Frame2*DropDownList1" ''( 7))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame4*Frame2*DropDownList1")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame4*Frame2*RealEntry3" ''( %s))',num2str(rhop));
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame5*Frame2*RealEntry3" ''( %s))',num2str(cpp));
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame6*Frame2*DropDownList1" ''( 5))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame6*Frame2*DropDownList1")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "User-Defined Functions*Frame1*List1" ''( 0))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "User-Defined Functions*Frame1*List1")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "User-Defined Functions*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Question*OK")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame7*Frame2*DropDownList1" ''( 11))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame7*Frame2*DropDownList1")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "User-Defined Functions*Frame1*List1" ''( 1))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "User-Defined Functions*Frame1*List1")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "User-Defined Functions*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame21*Frame2*RealEntry3" ''( %s))',num2str(betap));
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Create/Edit Materials*PanelButtons*PushButton1(Change/Create)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Create/Edit Materials*PanelButtons*PushButton1(Close)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Materials*Frame1*Table1*Frame1*List1(Materials)" ''( 4))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Materials*Frame1*Table1*Frame1*List1(Materials)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Materials*Frame1*Table1*ButtonBox2*PushButton1(Create/Edit)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Create/Edit Materials*Frame2(Properties)*Table2(Properties)*Frame6*Frame2*RealEntry3" ''( %s))',num2str(k_potop(a5)*kp)); %Porous region thermal conductivity
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Create/Edit Materials*PanelButtons*PushButton1(Change/Create)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Create/Edit Materials*PanelButtons*PushButton1(Close)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "MenuBar*DefineMenu*Phases...")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Phases*Frame1*Table1*Frame1*List1(Phases)" ''( 0))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Phases*Frame1*Table1*Frame1*List1(Phases)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Phases*Frame1*Table1*ButtonBox2*PushButton1(Edit)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "phase-domain-2*Frame2*Table2*Frame1*Table1*DropDownList1(Phase Material)" ''( 0))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "phase-domain-2*Frame2*Table2*Frame1*Table1*DropDownList1(Phase Material)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "phase-domain-2*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Phases*Frame1*Table1*Frame1*List1(Phases)" ''( 1))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Phases*Frame1*Table1*Frame1*List1(Phases)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Phases*Frame1*Table1*ButtonBox2*PushButton1(Edit)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "phase-domain-3*Frame3(Properties)*Table3(Properties)*Frame1*Frame2*RealEntry3" ''( %se-009))',num2str(dp));
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "phase-domain-3*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Phases*Frame1*Table1*ButtonBox2*PushButton2(Interaction)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "interaction-domain-4*PanelButtons*PushButton2(Cancel)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Phases*PanelButtons*PushButton1(Close)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Setup|Cell Zone Conditions"))');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Boundary Conditions*Frame1*Table1*Frame1*List1(Zone)" ''( 0))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Boundary Conditions*Frame1*Table1*Frame1*List1(Zone)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Boundary Conditions*Frame1*Table1*Frame2*Table2*Frame4*Table4*ButtonBox1*PushButton1(Edit)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-toggle-button "fluid-3-1*Frame4*Table4*CheckButton3(Porous Zone)" #f)');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "fluid-3-1*Frame4*Table4*CheckButton3(Porous Zone)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "fluid-3-1*Frame5*Table5*Frame1*Frame5(Porous Zone)*Frame1*Table1*Frame1*Table1*Frame7(Fluid Porosity)*Table7(Fluid Porosity)*Frame1*Table1*RealEntry2(Porosity)" ''( %s))',num2str(eps));
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "fluid-3-1*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Boundary Conditions*Frame1*Table1*Frame2*Table2*DropDownList1(Phase)" ''( 1))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Boundary Conditions*Frame1*Table1*Frame2*Table2*DropDownList1(Phase)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Boundary Conditions*Frame1*Table1*Frame1*List1(Zone)" ''( 0))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Boundary Conditions*Frame1*Table1*Frame1*List1(Zone)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Boundary Conditions*Frame1*Table1*Frame2*Table2*Frame4*Table4*ButtonBox1*PushButton1(Edit)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "fluid-3-2*Frame5*Table5*Frame1*Frame5(Porous Zone)*Frame1*Table1*Frame1*Table1*Frame4(Viscous Resistance (Inverse Absolute Permeability))*Table4(Viscous Resistance (Inverse Absolute Permeability))*Frame1*Table1*RealEntry2(Direction-1)" ''( %s))',num2str(viscous_resistance(a1,a2)));
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "fluid-3-2*Frame5*Table5*Frame1*Frame5(Porous Zone)*Frame1*Table1*Frame1*Table1*Frame4(Viscous Resistance (Inverse Absolute Permeability))*Table4(Viscous Resistance (Inverse Absolute Permeability))*Frame2*Table2*RealEntry2(Direction-2)" ''( %s))',num2str(viscous_resistance(a1,a2)));
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "fluid-3-2*Frame5*Table5*Frame1*Frame5(Porous Zone)*Frame1*Table1*Frame1*Table1*Frame5(Inertial Resistance)*Table5(Inertial Resistance)*Frame2*Table2*RealEntry2(Direction-1)" ''( %s))',num2str(c2(a1,a2)));
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "fluid-3-2*Frame5*Table5*Frame1*Frame5(Porous Zone)*Frame1*Table1*Frame1*Table1*Frame5(Inertial Resistance)*Table5(Inertial Resistance)*Frame3*Table3*RealEntry2(Direction-2)" ''( %s))',num2str(c2(a1,a2)));
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "fluid-3-2*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Boundary Conditions*Frame1*Table1*Frame2*Table2*DropDownList1(Phase)" ''( 2))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Boundary Conditions*Frame1*Table1*Frame2*Table2*DropDownList1(Phase)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Boundary Conditions*Frame1*Table1*Frame1*List1(Zone)" ''( 0))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Boundary Conditions*Frame1*Table1*Frame1*List1(Zone)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Boundary Conditions*Frame1*Table1*Frame2*Table2*Frame4*Table4*ButtonBox1*PushButton1(Edit)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "fluid-3-3*Frame5*Table5*Frame1*Frame5(Porous Zone)*Frame1*Table1*Frame1*Table1*Frame4(Viscous Resistance (Inverse Absolute Permeability))*Table4(Viscous Resistance (Inverse Absolute Permeability))*Frame1*Table1*RealEntry2(Direction-1)" ''( %s))',num2str(viscous_resistance(a1,a2)));
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "fluid-3-3*Frame5*Table5*Frame1*Frame5(Porous Zone)*Frame1*Table1*Frame1*Table1*Frame4(Viscous Resistance (Inverse Absolute Permeability))*Table4(Viscous Resistance (Inverse Absolute Permeability))*Frame2*Table2*RealEntry2(Direction-2)" ''( %s))',num2str(viscous_resistance(a1,a2)));
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "fluid-3-3*Frame5*Table5*Frame1*Frame5(Porous Zone)*Frame1*Table1*Frame1*Table1*Frame5(Inertial Resistance)*Table5(Inertial Resistance)*Frame2*Table2*RealEntry2(Direction-1)" ''( %s))',num2str(c2(a1,a2)));
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "fluid-3-3*Frame5*Table5*Frame1*Frame5(Porous Zone)*Frame1*Table1*Frame1*Table1*Frame5(Inertial Resistance)*Table5(Inertial Resistance)*Frame3*Table3*RealEntry2(Direction-2)" ''( %s))',num2str(c2(a1,a2)));
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "fluid-3-3*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Boundary Conditions*Frame1*Table1*Frame2*Table2*DropDownList1(Phase)" ''( 0))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Boundary Conditions*Frame1*Table1*Frame2*Table2*DropDownList1(Phase)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-toggle-button "Boundary Conditions*Frame1*Table1*ButtonBox3(Porous Formulation)*Physical Velocity" #f)');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Boundary Conditions*Frame1*Table1*ButtonBox3(Porous Formulation)*Physical Velocity")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Boundary Conditions*Frame1*Table1*Frame2*Table2*Frame4*Table4*Frame2*Table2*PushButton2(Operating Conditions)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Operating Conditions*Frame1*Frame3(Gravity)*Table3(Gravity)*Frame3(Boussinesq Parameters)*Table3(Boussinesq Parameters)*RealEntry1(Operating Temperature)" ''( %s))',num2str((Th+Tc)/2));
                    fprintf(jou,'\r\n(cx-gui-do cx-set-toggle-button "Operating Conditions*Frame1*Frame3(Gravity)*Table3(Gravity)*Frame4(Variable-Density Parameters)*Table4(Variable-Density Parameters)*ToggleBox1*CheckButton1(Specified Operating Density)" #f)');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Operating Conditions*Frame1*Frame3(Gravity)*Table3(Gravity)*Frame4(Variable-Density Parameters)*Table4(Variable-Density Parameters)*ToggleBox1*CheckButton1(Specified Operating Density)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Operating Conditions*Frame1*Frame3(Gravity)*Table3(Gravity)*Frame4(Variable-Density Parameters)*Table4(Variable-Density Parameters)*RealEntry3(Operating Density)" ''( 995))');     %Operating density-if needed add str2num
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Operating Conditions*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Setup|Boundary Conditions"))');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Boundary Conditions*Frame1*Table1*Frame1*List1(Zone)" ''( 4))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Boundary Conditions*Frame1*Table1*Frame1*List1(Zone)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Boundary Conditions*Frame1*Table1*Frame2*Table2*Frame4*Table4*ButtonBox1*PushButton1(Edit)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-toggle-button "wall-5-1*Frame4*Frame3(Thermal)*Frame1*Frame1(Thermal Conditions)*Frame1*ToggleBox1*Temperature" #f)');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "wall-5-1*Frame4*Frame3(Thermal)*Frame1*Frame1(Thermal Conditions)*Frame1*ToggleBox1*Temperature")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "wall-5-1*Frame4*Frame3(Thermal)*Frame1*Frame1(Thermal Conditions)*Frame7*Table7*Frame1*Table1*RealEntry2(Temperature)" ''( %s))',num2str(Th));
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "wall-5-1*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Boundary Conditions*Frame1*Table1*Frame1*List1(Zone)" ''( 3))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Boundary Conditions*Frame1*Table1*Frame1*List1(Zone)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Boundary Conditions*Frame1*Table1*Frame2*Table2*Frame4*Table4*ButtonBox1*PushButton1(Edit)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-toggle-button "wall-4-1*Frame4*Frame3(Thermal)*Frame1*Frame1(Thermal Conditions)*Frame1*ToggleBox1*Temperature" #f)');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "wall-4-1*Frame4*Frame3(Thermal)*Frame1*Frame1(Thermal Conditions)*Frame1*ToggleBox1*Temperature")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "wall-4-1*Frame4*Frame3(Thermal)*Frame1*Frame1(Thermal Conditions)*Frame7*Table7*Frame1*Table1*RealEntry2(Temperature)" ''( %s))',num2str(Tc));
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "wall-4-1*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Solution|Solution Methods"))');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Solution|Solution Controls"))');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Solution|Monitors"))');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Monitors*Frame1*Table1*Frame1*List1(Residuals, Statistic and Force Monitors)" ''( 0))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Monitors*Frame1*Table1*Frame1*List1(Residuals, Statistic and Force Monitors)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Monitors*Frame1*Table1*Frame2*Table2*PushButton2(Edit)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Residual Monitors*Frame1*Table1*Frame2*Table2*Frame1(Equations)*Table1(Equations)*RealEntry11" ''( 0.0001))');   %continuity residuals
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Residual Monitors*Frame1*Table1*Frame2*Table2*Frame1(Equations)*Table1(Equations)*RealEntry17" ''( 5e-005))');   %velocity residuals
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Residual Monitors*Frame1*Table1*Frame2*Table2*Frame1(Equations)*Table1(Equations)*RealEntry23" ''( 5e-005))');   %velocity residuals
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Residual Monitors*Frame1*Table1*Frame2*Table2*Frame1(Equations)*Table1(Equations)*RealEntry29" ''( 1e-009))');   %Energy residuals
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Residual Monitors*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Monitors*Frame1*Table1*Frame4*Table4*PushButton1(Create)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Surface Monitor*Frame1*Table1*Frame2*Table2*DropDownList3(Field Variable)" ''( 5))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Surface Monitor*Frame1*Table1*Frame2*Table2*DropDownList3(Field Variable)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Surface Monitor*Frame1*Table1*Frame2*Table2*DropDownList4" ''( 4))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Surface Monitor*Frame1*Table1*Frame2*Table2*DropDownList4")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Surface Monitor*Frame1*Table1*Frame2*Table2*Frame6*Table6*Frame1*List1(Surfaces)" ''( 4))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Surface Monitor*Frame1*Table1*Frame2*Table2*Frame6*Table6*Frame1*List1(Surfaces)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-toggle-button "Surface Monitor*Frame1*Table1*Frame1*Table1*Frame2(Options)*Table2(Options)*CheckButton2(Plot)" #f)');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Surface Monitor*Frame1*Table1*Frame1*Table1*Frame2(Options)*Table2(Options)*CheckButton2(Plot)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Surface Monitor*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n/solve/set/expert yes no no no yes no ');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Solution|Solution Initialization"))');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Solution Initialization*Frame1*Table1*DropDownList1(Compute from)" ''( 0))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Solution Initialization*Frame1*Table1*DropDownList1(Compute from)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-toggle-button "Solution Initialization*Frame1*Table1*Frame2(Reference Frame)*ToggleBox2(Reference Frame)*Absolute" #f)');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Solution Initialization*Frame1*Table1*Frame2(Reference Frame)*ToggleBox2(Reference Frame)*Absolute")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Solution Initialization*Frame1*Table1*Frame6(Initial Values)*Table6(Initial Values)*RealEntry5(phase-2 Volume Fraction)" ''( %s))',num2str(phi));
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Solution Initialization*Frame1*Table1*Frame6(Initial Values)*Table6(Initial Values)*RealEntry5(phase-2 Volume Fraction)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Solution Initialization*Frame1*Table1*ButtonBox8*PushButton1(Initialize)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Solution|Run Calculation"))');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-integer-entry "Run Calculation*Frame1*Table1*IntegerEntry9(Number of Iterations)" 500)'); %First iteration number
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Run Calculation*Frame1*Table1*IntegerEntry9(Number of Iterations)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Run Calculation*Frame1*Table1*PushButton21(Calculate)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "GraphicsArea*GraphicsView1*DropDownList1" ''( 0))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "GraphicsArea*GraphicsView1*DropDownList1")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Information*OK")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Solution|Monitors"))');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Monitors*Frame1*Table1*Frame1*List1(Residuals, Statistic and Force Monitors)" ''( 0))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Monitors*Frame1*Table1*Frame1*List1(Residuals, Statistic and Force Monitors)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Monitors*Frame1*Table1*Frame2*Table2*PushButton2(Edit)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Residual Monitors*Frame1*Table1*Frame2*Table2*Frame1(Equations)*Table1(Equations)*RealEntry17" ''( 1e-006))'); %second velocity residuals
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Residual Monitors*Frame1*Table1*Frame2*Table2*Frame1(Equations)*Table1(Equations)*RealEntry23" ''( 1e-006))'); %second velocity residuals
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Residual Monitors*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Solution|Solution Controls"))');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Solution Controls*Frame1*Table1*Frame5(Under-Relaxation Factors)*Table5(Under-Relaxation Factors)*RealEntry2(Density)" ''( 0.9))'); %Density relaxation factor
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Solution Controls*Frame1*Table1*Frame5(Under-Relaxation Factors)*Table5(Under-Relaxation Factors)*RealEntry2(Density)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Solution Controls*Frame1*Table1*Frame5(Under-Relaxation Factors)*Table5(Under-Relaxation Factors)*RealEntry3(Body Forces)" ''( 0.9))'); %Body force relaxation factor
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Solution Controls*Frame1*Table1*Frame5(Under-Relaxation Factors)*Table5(Under-Relaxation Factors)*RealEntry3(Body Forces)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Solution Controls*Frame1*Table1*Frame5(Under-Relaxation Factors)*Table5(Under-Relaxation Factors)*RealEntry4(Momentum)" ''( 0.4))'); %Momentum relaxation factor
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Solution Controls*Frame1*Table1*Frame5(Under-Relaxation Factors)*Table5(Under-Relaxation Factors)*RealEntry4(Momentum)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Solution|Run Calculation"))');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Solution Controls*Frame1*Table1*Frame5(Under-Relaxation Factors)*Table5(Under-Relaxation Factors)*RealEntry7(Energy)" ''( 0.9))'); %Energy relaxation factor
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Solution Controls*Frame1*Table1*Frame5(Under-Relaxation Factors)*Table5(Under-Relaxation Factors)*RealEntry7(Energy)")');

                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "MenuBar*User-DefinedSubMenu*Memory...")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-integer-entry "User-Defined Memory*Frame1*Table1*IntegerEntry2" 1)');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "User-Defined Memory*Frame1*Table1*IntegerEntry2")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "User-Defined Memory*Frame1*Table1*IntegerEntry2")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-integer-entry "User-Defined Memory*Frame1*Table1*IntegerEntry2" 2)');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "User-Defined Memory*Frame1*Table1*IntegerEntry2")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "User-Defined Memory*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "MenuBar*User-DefinedSubMenu*Execute on Demand...")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-selections "Execute on Demand*DropDownList1(Execute on Demand)" ''( 1))');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Execute on Demand*DropDownList1(Execute on Demand)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Execute on Demand*PanelButtons*PushButton1(OK)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Execute on Demand*PanelButtons*PushButton2(Cancel)")'); 
                    fprintf(jou,'\r\n(cx-gui-do cx-set-integer-entry "Run Calculation*Frame1*Table1*IntegerEntry9(Number of Iterations)" 1500)'); %Number of second iteration
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Run Calculation*Frame1*Table1*IntegerEntry9(Number of Iterations)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Run Calculation*Frame1*Table1*PushButton21(Calculate)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Information*OK")'); 
                    
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Solution|Solution Controls"))');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Solution Controls*Frame1*Table1*Frame5(Under-Relaxation Factors)*Table5(Under-Relaxation Factors)*RealEntry2(Density)" ''( 0.85))'); %Density relaxation factor
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Solution Controls*Frame1*Table1*Frame5(Under-Relaxation Factors)*Table5(Under-Relaxation Factors)*RealEntry2(Density)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Solution Controls*Frame1*Table1*Frame5(Under-Relaxation Factors)*Table5(Under-Relaxation Factors)*RealEntry3(Body Forces)" ''( 0.87))'); %Body force relaxation factor
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Solution Controls*Frame1*Table1*Frame5(Under-Relaxation Factors)*Table5(Under-Relaxation Factors)*RealEntry3(Body Forces)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Solution Controls*Frame1*Table1*Frame5(Under-Relaxation Factors)*Table5(Under-Relaxation Factors)*RealEntry4(Momentum)" ''( 0.1))'); %Momentum relaxation factor
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Solution Controls*Frame1*Table1*Frame5(Under-Relaxation Factors)*Table5(Under-Relaxation Factors)*RealEntry4(Momentum)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-list-tree-selections "NavigationPane*List_Tree1" (list "Solution|Run Calculation"))');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-real-entry-list "Solution Controls*Frame1*Table1*Frame5(Under-Relaxation Factors)*Table5(Under-Relaxation Factors)*RealEntry7(Energy)" ''( 0.85))'); %Energy relaxation factor
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Solution Controls*Frame1*Table1*Frame5(Under-Relaxation Factors)*Table5(Under-Relaxation Factors)*RealEntry7(Energy)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-integer-entry "Run Calculation*Frame1*Table1*IntegerEntry9(Number of Iterations)" 1000)'); %Number of third iteration
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Run Calculation*Frame1*Table1*IntegerEntry9(Number of Iterations)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Run Calculation*Frame1*Table1*PushButton21(Calculate)")');
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Information*OK")'); 
                    
                    fprintf(jou,'\r\n/report/surface-integrals/integral t_h () mixture heat-flux y "C:\\Users\\---\\finalruns\\textresult"'); %Location of exporting results
                    fprintf(jou,'\r\n/report/volume-integrals/volume-integral inner outer () mixture udm-0 y "C:\\Users\\---\\finalruns\\sgen0"'); %Location of exporting entropy gen 0
                    fprintf(jou,'\r\n/report/volume-integrals/volume-integral inner outer () mixture udm-1 y "C:\\Users\\---\\finalruns\\sgen1"'); %Location of exporting entropy gen 1
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "MenuBar*WriteSubMenu*Case & Data...")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-text-entry "Select File*FilterText" "C:\\Users\\---\\finalruns\\*")'); %Location of case and data
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Select File*Apply")');
                    fprintf(jou,'\r\n(cx-gui-do cx-set-text-entry "Select File*Text" "%s%s%s%s%s")',num2str(a1),num2str(a2),num2str(a3),num2str(a4),num2str(a5)); %Name of case&data files: 5digit number that sequently shows Ra,Da,n,partial filled,ration of k matrix indicies
                    fprintf(jou,'\r\n(cx-gui-do cx-activate-item "Select File*OK")');
                    fprintf(jou,'\r\nexit');                    
                    fclose(jou);
                    %% Run Fluent with journal
                    !"C:\Program Files\ANSYS Inc\---\fluent.exe"  2ddp -wait -i "C:\---\j.jou"
                    %% Saving results
                    result(a1,a2,a3,a4,a5)= reportdigit('C:\Users\---\textresult'); %Total wall heat flux: getting the number from text report
                    Nusselt(a1,a2,a3,a4,a5)=result(a1,a2,a3,a4,a5)/(kf*(Th-Tc));
                    
                    s0(a1,a2,a3,a4,a5)= reportdigit('C:\Users\---\sgen0'); %Volume intergral for entropy generation: getting the number from text report
                    s1(a1,a2,a3,a4,a5)= reportdigit('C:\Users\---\sgen1'); %Volume intergral for entropy generation: getting the number from text report
                    stot=s0+s1;
                    
                    delete('C:\Users\---\textresult');
                    delete('C:\Users\---\sgen0');
                    delete('C:\Users\---\sgen1');
                end
            end
        end
    end
end
save('variables234','result','Nusselt','s0','s1','stot');
