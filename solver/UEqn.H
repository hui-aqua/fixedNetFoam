// Solve the Momentum equation
const scalar porosity(Nettings.porosity());

MRF.correctBoundaryVelocity(U);


tmp<fvVectorMatrix> tUEqn
(
    1/porosity*fvm::ddt(U) + fvm::div(phi, U)
  + 1/porosity*MRF.DDt(U)
  + 1/porosity*turbulence->divDevReff(U)
 ==
    fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

Nettings.addResistance(U,Fh,UEqn,mesh);


fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}

// Nettings.updateVelocity(mesh,UEqn);
