

float Ele_Eff_Area(float SCeta_stuff){

    if(fabs(SCeta_stuff) >= 0.0    && fabs(SCeta_stuff) < 1.    ){return 0.1752;}
    if(fabs(SCeta_stuff) >= 1.     && fabs(SCeta_stuff) < 1.479 ){return 0.1862;}
    if(fabs(SCeta_stuff) >= 1.479  && fabs(SCeta_stuff) < 2.0   ){return 0.1411;}
    if(fabs(SCeta_stuff) >= 2.0    && fabs(SCeta_stuff) < 2.2   ){return 0.1534;}
    if(fabs(SCeta_stuff) >= 2.2    && fabs(SCeta_stuff) < 2.3   ){return 0.1903;}
    if(fabs(SCeta_stuff) >= 2.3    && fabs(SCeta_stuff) < 2.4   ){return 0.2243;}
    if(fabs(SCeta_stuff) >= 2.4                                 ){return 0.2687;}

    else{
        std::cout << "Did not get Eff Area value" << std::endl;
        return 0.0;
    } 
    
}
