void main()
{

    double e_tank_lox_1 = 4.700519093954655e-06;
    double D_tank_lox_1 = 4.73994357405418;
    double L_tank_lox_1 = 5.191133141978805;
    double L_total_lox_1 = 9.931076716032987;


    double e_tank_h2_1 = 0.0036992490696796704;
    double D_tank_h2_1 = 4.73994357405418;
    double L_tank_h2_1 = 19.320043264457187;
    double L_total_h2_1 = 24.059986838511367;


    double e_stage_1 = 0.00310668666905349;
    double L_stage_1 = 35.082461518499386;
    double D_stage_1 = 4.770137703943644;


    double e_nozzle_1 = 231.66123468830767;
    double h_nozzle_1 = 3.3408356755555726;
    double D_throat_1 = 0.4921995827153984;
    double D_exit_1 = 2.282548025335798;
    double e_cc_1 = 231.66123468830767;
    double h_cone_1 = 0.10387575812268772;
    double D_cc_1 = 0.6999510989607738;
    double h_cc_1 = 0.4004919283043367;



/////////////////////////////////////////////////////////////
//////// LOX tank //////////////////////

    string tank_lox_1_id = AddGeom( "STACK" );

    CutXSec( tank_lox_1_id, 1 );

    string tank_lox_1_xsurf_id = GetXSecSurf( tank_lox_1_id, 0 );



    SetParmVal(tank_lox_1_id, "XDelta", "XSec_1", D_tank_lox_1 / 2);
    SetParmVal(tank_lox_1_id, "XDelta", "XSec_2", L_tank_lox_1); 
    SetParmVal(tank_lox_1_id, "XDelta", "XSec_3",D_tank_lox_1 / 2); 

    SetParmVal(tank_lox_1_id, "X_Rel_Location", "XForm", -(L_total_lox_1 -0.5));
    SetParmVal(tank_lox_1_id, "Y_Rel_Location", "XForm", 0);


    ChangeXSecShape( tank_lox_1_xsurf_id, 1, XS_CIRCLE);    
    ChangeXSecShape( tank_lox_1_xsurf_id, 2, XS_CIRCLE);
 

    string tank_lox_1_xsec_id0 = GetXSec( tank_lox_1_xsurf_id, 0 );
    string tank_lox_1_xsec_id1 = GetXSec( tank_lox_1_xsurf_id, 1 );
    string tank_lox_1_xsec_id2 = GetXSec( tank_lox_1_xsurf_id, 2 );
    string tank_lox_1_xsec_id3 = GetXSec( tank_lox_1_xsurf_id, 3 );
 


    SetParmVal( GetXSecParm( tank_lox_1_xsec_id1, "Circle_Diameter"), D_tank_lox_1 );
    SetParmVal( GetXSecParm( tank_lox_1_xsec_id2, "Circle_Diameter"),D_tank_lox_1);



    SetParmVal( GetParm(tank_lox_1_id, "TopLAngle", "XSec_3" ),-90.0 );
    SetParmVal( GetParm(tank_lox_1_id, "RightLAngle", "XSec_3" ),-90.0 );

    SetParmVal( GetParm(tank_lox_1_id, "TopLAngle", "XSec_0" ), 90.0 );
    SetParmVal( GetParm(tank_lox_1_id, "RightLAngle", "XSec_0" ), 90.0 );


/////////////////////////////////////////////////////////////////////////
///////////////          H2 tank        /////////////////////////////////

    string tank_h2_1_id = AddGeom( "STACK" );

    CutXSec( tank_h2_1_id, 1 );

    string tank_h2_1_xsurf_id = GetXSecSurf( tank_h2_1_id, 0 );



    SetParmVal(tank_h2_1_id, "XDelta", "XSec_1", D_tank_h2_1 / 2);
    SetParmVal(tank_h2_1_id, "XDelta", "XSec_2", L_tank_h2_1); 
    SetParmVal(tank_h2_1_id, "XDelta", "XSec_3",D_tank_h2_1 / 2); 

    SetParmVal(tank_h2_1_id, "X_Rel_Location", "XForm", -(L_stage_1-0.5));
    SetParmVal(tank_h2_1_id, "Y_Rel_Location", "XForm", 0);


    ChangeXSecShape( tank_h2_1_xsurf_id, 1, XS_CIRCLE);    
    ChangeXSecShape( tank_h2_1_xsurf_id, 2, XS_CIRCLE);
 

    string tank_h2_1_xsec_id0 = GetXSec(  tank_h2_1_xsurf_id, 0 );
    string tank_h2_1_xsec_id1 = GetXSec(  tank_h2_1_xsurf_id, 1 );
    string tank_h2_1_xsec_id2 = GetXSec(  tank_h2_1_xsurf_id, 2 );
    string tank_h2_1_xsec_id3 = GetXSec(  tank_h2_1_xsurf_id, 3 );
 


    SetParmVal( GetXSecParm( tank_h2_1_xsec_id1, "Circle_Diameter"), D_tank_h2_1 );
    SetParmVal( GetXSecParm( tank_h2_1_xsec_id2, "Circle_Diameter"), D_tank_h2_1);



    SetParmVal( GetParm(tank_h2_1_id, "TopLAngle", "XSec_3" ),-90.0 );
    SetParmVal( GetParm(tank_h2_1_id, "RightLAngle", "XSec_3" ),-90.0 );

    SetParmVal( GetParm(tank_h2_1_id, "TopLAngle", "XSec_0" ), 90.0 );
    SetParmVal( GetParm(tank_h2_1_id, "RightLAngle", "XSec_0" ), 90.0 );

///////////////////////////////////////////////////////////////////////
/////////////////////////  COMBUSTION CHAMBER //////////////////////////

    string cc_1_id = AddGeom( "STACK" );

    CutXSec( cc_1_id, 1 );

    string cc_1_xsurf_id = GetXSecSurf( cc_1_id, 0 );



    SetParmVal(cc_1_id, "XDelta", "XSec_1", 0);
    SetParmVal(cc_1_id, "XDelta", "XSec_2", h_cc_1); 
    SetParmVal(cc_1_id, "XDelta", "XSec_3",h_cone_1); 


    SetParmVal(cc_1_id, "X_Rel_Location", "XForm", 1);
    SetParmVal(cc_1_id, "Y_Rel_Location", "XForm", 0);


    ChangeXSecShape( cc_1_xsurf_id, 1, XS_CIRCLE);    
    ChangeXSecShape( cc_1_xsurf_id, 2, XS_CIRCLE);
    ChangeXSecShape( cc_1_xsurf_id, 3, XS_CIRCLE);


    string cc_1_xsec_id0 = GetXSec( cc_1_xsurf_id, 0 );
    string cc_1_xsec_id1 = GetXSec( cc_1_xsurf_id, 1 );
    string cc_1_xsec_id2 = GetXSec( cc_1_xsurf_id, 2 );
    string cc_1_xsec_id3 = GetXSec( cc_1_xsurf_id, 3 );




    SetParmVal( GetXSecParm( cc_1_xsec_id1, "Circle_Diameter"), D_cc_1 );
    SetParmVal( GetXSecParm( cc_1_xsec_id2, "Circle_Diameter"),D_cc_1);
    SetParmVal( GetXSecParm( cc_1_xsec_id3, "Circle_Diameter"),D_throat_1);

//////////////////////////////////////////////////////////////////
////////////////////////// NOZZLE  /////////////////////////////////

    string nozzle_1_id = AddGeom( "STACK" );

    CutXSec( nozzle_1_id, 1 );

    string nozzle_1_xsurf_id = GetXSecSurf( nozzle_1_id,0 );



    SetParmVal(nozzle_1_id, "XDelta", "XSec_0", 0);
    SetParmVal(nozzle_1_id, "XDelta", "XSec_1",h_nozzle_1);
    SetParmVal(nozzle_1_id, "XDelta", "XSec_2",0);  

    CutXSec( nozzle_1_id, 3 );

    SetParmVal(nozzle_1_id, "X_Rel_Location", "XForm", h_cc_1 + h_cone_1+1);
    SetParmVal(nozzle_1_id, "Y_Rel_Location", "XForm", 0);


    ChangeXSecShape( nozzle_1_xsurf_id, 0, XS_CIRCLE);    
    ChangeXSecShape( nozzle_1_xsurf_id, 1, XS_CIRCLE);
    ChangeXSecShape( nozzle_1_xsurf_id, 2, XS_CIRCLE);



    string nozzle_1_xsec_id1 = GetXSec( nozzle_1_xsurf_id, 0 );
    string nozzle_1_xsec_id2 = GetXSec( nozzle_1_xsurf_id, 1 );
    string nozzle_1_xsec_id3 = GetXSec( nozzle_1_xsurf_id, 2 );
    string nozzle_1_xsec_id4 = GetXSec( nozzle_1_xsurf_id, 3 );




    SetParmVal( GetXSecParm( nozzle_1_xsec_id1, "Circle_Diameter"), D_throat_1 );
    SetParmVal( GetXSecParm( nozzle_1_xsec_id2, "Circle_Diameter"),D_exit_1);
    SetParmVal( GetXSecParm( nozzle_1_xsec_id3, "Circle_Diameter"),D_exit_1);


///////////////////////////////////////////////////////////////////////////
//////////////////////// FIRST STAGE //////////////////////////////////////



   string stage_1_id = AddGeom( "STACK" );

    CutXSec( stage_1_id, 1 );

    string stage_1_xsurf_id = GetXSecSurf( stage_1_id, 0 );



    SetParmVal(stage_1_id, "XDelta", "XSec_1", 0.5);
    SetParmVal(stage_1_id, "XDelta", "XSec_2", L_stage_1); 
    SetParmVal(stage_1_id, "XDelta", "XSec_3",0.5); 

    SetParmVal(stage_1_id, "X_Rel_Location", "XForm", -L_stage_1);
    SetParmVal(stage_1_id, "Y_Rel_Location", "XForm", 0);


    ChangeXSecShape( stage_1_xsurf_id, 1, XS_CIRCLE);    
    ChangeXSecShape( stage_1_xsurf_id, 2, XS_CIRCLE);
 

    string stage_1_xsec_id0 = GetXSec(  stage_1_xsurf_id, 0 );
    string stage_1_xsec_id1 = GetXSec(  stage_1_xsurf_id, 1 );
    string stage_1_xsec_id2 = GetXSec(  stage_1_xsurf_id, 2 );
    string stage_1_xsec_id3 = GetXSec(  stage_1_xsurf_id, 3 );
 


    SetParmVal( GetXSecParm( stage_1_xsec_id1, "Circle_Diameter"), D_stage_1);
    SetParmVal( GetXSecParm( stage_1_xsec_id2, "Circle_Diameter"), D_stage_1);



    SetParmVal( GetParm(stage_1_id, "TopLAngle", "XSec_3" ),-90.0 );
    SetParmVal( GetParm(stage_1_id, "RightLAngle", "XSec_3" ),-90.0 );

    SetParmVal( GetParm(stage_1_id, "TopLAngle", "XSec_0" ), 90.0 );
    SetParmVal( GetParm(stage_1_id, "RightLAngle", "XSec_0" ), 90.0 );




/////////////////////////////////////
/////////////////////////////////////

 
    WriteVSPFile( "LAST.vsp3", SET_ALL ); // Write To File


}