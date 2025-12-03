

const PROJECT_ROOT = "C:/Users/bettif/Documents/GitHub/abaqus-substructuring";

function main() {
    console.log("Starting the piecewise solver...");

    console.log("Opening the model...");
    var model = Spck.openModel(PROJECT_ROOT + "/Wheelset.spck");

    model.createBody("$B_FlexTrack", false);
    var obj2 = model.findElement("$B_FlexTrack");
    obj2.flx.file.src = PROJECT_ROOT + "/SplineSegment01_FREQ1.fbi";
    obj2.type.src = 9;
    obj2.createPrim("$P_FlexTrack_Flex");
    var obj3 = model.findElement("$P_FlexTrack_Flex");
    obj3.par(1).src = 0.11140086376720559;
    obj3.par(2).src = 0.14853448502294081;
    obj3.par(3).src = 0.18566810627867597;
    obj3.type.src = 28;
    obj3.par(2).src = 1;
    obj3.par(3).src = 1;
    obj3.mp.incl.src = 0;

    obj2.createMarker("$M_FlexTrack_fe_reference");
    var obj4 = model.findElement("$M_FlexTrack_fe_reference");
    obj4.flx.type.src = 4;
    obj4.pos(1).src = -0.01;
    obj4.pos(2).src = 0.32000000000000001;
    obj2.flx.ref.src = "$M_FlexTrack_fe_reference";
    
    model.destroy(model.findElement("$P_FlexTrack"));

    model.createForce("$F_StopIntegration");
    var obj6 = model.findElement("$F_StopIntegration");
    obj6.type.src = 234;
    obj6.par(1).src = 50;
    obj6.par(3).src = "$J_Wheelset";
    obj6.par(4).src = 1;

    var obj = model.findElement("$SLV_SolverSettings");
    obj.integ.tend.time.src = 1000;

    for (var i = 1; i <= 3; i++) {
        num_segment = i.toString().padStart(2, '0');
        obj2.flx.file.src = PROJECT_ROOT + "/SplineSegment" + num_segment + "_FREQ1.fbi";

        for (let k = 0; k < 1000; k++) {
            obj2.flx.st.active(k).src = 1;
            obj2.flx.st.prestress(k).src = 0;
            obj2.flx.st.scal.stiff(k).src = 1;
            obj2.flx.st.dep(k).src = 2;
            obj2.flx.st.vel(k).src = 0;
            obj2.flx.st.pos(k).src = 0;
        }

        obj6.par(1).src = 50 * i;
        if (i > 1) {
            obj6.par(1).src = obj6.par(1).src - 20;
        }
        if (i == 1) {
            Spck.Slv.preld(model);
            Spck.currentModel.vehicle.startvel.src = "5 m/s";
            Spck.Slv.vehicleVelocities(model);
        }
        Spck.Slv.integMeas(model);
    }
}