/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
 
 No.

\*---------------------------------------------------------------------------*/

#include "netPanel.H"
#include "volFields.H"
#include "SortableList.H"
// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::netPanel::calcNorm(
        const point &pointI,
        const point &pointII,
        const point &pointIII) const
{
    const vector a(pointI - pointII);
    const vector b(pointI - pointIII);
    const vector norm(a ^ b / (mag(a ^ b) + SMALL));
    return norm; 
}

Foam::scalar Foam::netPanel::calcTheta(
        const point &pointI,
        const point &pointII,
        const point &pointIII,
        const vector &fluidVelocity) const
{
    const vector a(pointI - pointII);
    const vector b(pointI - pointIII);
    const scalar theta((a&b)/(mag(a)*mag(b)));  

    return theta;  
}


Foam::scalar Foam::netPanel::calcArea(
        const point &pointI,
        const point &pointII,
        const point &pointIII) const
{
    const vector a(pointI - pointII);
    const vector b(pointI - pointIII);
    return 0.5 * mag(a ^ b);
}


Foam::scalar Foam::netPanel::calcDist(
        const point &pointI,
        const point &pointII) const
{
    return mag(pointI - pointII);
}


Foam::scalar Foam::netPanel::calcDistanceFromPoint2Panel(
    const point &x,
    const vector &structuralElementi) const
{
    const point point_a = structuralPositions_memb[structuralElementi[0]];
    const point point_b = structuralPositions_memb[structuralElementi[1]];
    const point point_c = structuralPositions_memb[structuralElementi[2]];
    vector panelNorm = calcNorm(point_a, point_b, point_c); // a unit vector to indicate the normal
    scalar dis(mag((x - point_a) & panelNorm));
    //    Info << "The distance from point to net panel is "<< dis << " m." << endl;
    return dis;
}

bool Foam::netPanel::isInPorous_line(
    const point &x,
    const vector &point_a,
    const vector &point_b,
    const vector &point_c) const
{
    bool result(false);

    // Method to enhance the mesh line ()
    const Foam::scalar line1(mag(point_a - point_b));
    const Foam::scalar line2(mag(point_a - point_c));
    const Foam::scalar line3(mag(point_b - point_c));
    Foam::scalar d1(0);
    Foam::scalar d2(0);
    Foam::scalar dis2line(0);
    if (line1 > line2 and line1 > line3)
    {
        // based on line 2
        d1 = (mag(point_a - x));
        d2 = (mag(point_c - x));
        dis2line = mag(calcArea(point_a, point_c, x) * 2.0 / (line2 + Foam::SMALL));
        if (dis2line <= thickness_memb * 0.5 * ropeEnhance_memb and max(d1, d2) <= (thickness_memb * 0.5 * ropeEnhance_memb + line2))
        {
            result = true;
        }
        else
        {
            // based on line3
            d1 = (mag(point_b - x));
            d2 = (mag(point_c - x));
            dis2line = (mag(calcArea(point_b, point_c, x) * 2.0 / (line3 + Foam::SMALL)));
            if (dis2line <= thickness_memb * 0.5 * ropeEnhance_memb and max(d1, d2) <= (thickness_memb * 0.5 * ropeEnhance_memb + line3))
            {
                result = true;
            }
        }
    }
    else
    {
        if (line2 < line3)
        {
            //line 2
            d1 = (mag(point_a - x));
            d2 = (mag(point_c - x));
            dis2line = (mag(calcArea(point_a, point_c, x) * 2.0 / (line2 + Foam::SMALL)));
            if (dis2line <= thickness_memb * 0.5 * ropeEnhance_memb and max(d1, d2) <= (thickness_memb * 0.5 * ropeEnhance_memb + line2))
            {
                result = true;
            }
            else
            {
                //line 1
                d1 = (mag(point_a - x));
                d2 = (mag(point_b - x));
                dis2line = (mag(calcArea(point_a, point_b, x) * 2.0 / (line1 + Foam::SMALL)));
                if (dis2line <= thickness_memb * 0.5 * ropeEnhance_memb and max(d1, d2) <= (thickness_memb * 0.5 * ropeEnhance_memb + line1))
                {
                    result = true;
                }
            }
        }
        else
        {
            //line 3
            d1 = (mag(point_b - x));
            d2 = (mag(point_c - x));
            dis2line = (mag(calcArea(point_b, point_c, x) * 2.0 / (line3 + Foam::SMALL)));
            if (dis2line <= thickness_memb * 0.5 * ropeEnhance_memb and max(d1, d2) <= (thickness_memb * 0.5 * ropeEnhance_memb + line3))
            {
                result = true;
            }
            else
            {
                //line 1
                d1 = (mag(point_a - x));
                d2 = (mag(point_b - x));
                dis2line = (mag(calcArea(point_a, point_b, x) * 2.0 / (line1 + Foam::SMALL)));
                if (dis2line <= thickness_memb * 0.5 * ropeEnhance_memb and max(d1, d2) <= (thickness_memb * 0.5 * ropeEnhance_memb + line1))
                {
                    result = true;
                }
            }
        }
    }
    return result;
}

bool Foam::netPanel::isInPorousZone(
    const point &x,
    const vector &point_a,
    const vector &point_b,
    const vector &point_c) const
{
    bool result(false); // initial value

    vector panelNorm = calcNorm(point_a, point_b, point_c); // a unit vector to indicate the normal
    scalar dis(mag((x - point_a) & panelNorm));
    // define a const scalar as the distance between point x to net panel
    if (dis <= thickness_memb * 0.5) // distance is less than half thickness
    {
        scalar panelarea(calcArea(point_a, point_b, point_c));
        vector projectedPoint(0, 0, 0);        // initial the projected point is 0,0,0
        if (((x - point_a) & panelNorm) < 0.0) // on the side of normal vector
        {
            projectedPoint = (x + panelNorm * dis);
        }
        else
        {
            projectedPoint = (x - panelNorm * dis);
        }
        // projectedPiont is the projected point on the net panel
        scalar panelarea3(calcArea(point_a, point_b, projectedPoint) +
                          calcArea(point_a, projectedPoint, point_c) +
                          calcArea(projectedPoint, point_b, point_c)); //  the area of the three trigular shapes.
        if (panelarea3 <= Foam::SMALL + panelarea)
        {
            result = true;
        }
    }

    return result;
}




// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::netPanel::netPanel(
        const dictionary &netDict):
        // initial components
        netDict_memb(netDict),
        Sn_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("Sn"))),
        thickness_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("PorousMediaThickness"))),
        dw_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("twineDiameter"))),
        ropeEnhance_memb(readScalar(netDict_memb.subDict("NetInfo1").lookup("ropeEnhance")))     
        
{
    // creat the object.
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::netPanel::~netPanel()
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::netPanel::addResistance(
        const volVectorField &U,
        volVectorField &Fh,
        fvVectorMatrix &UEqn,
        // volScalarField &porosityField,
        const fvMesh &mesh) const
{
    Info << ">>> addResistance" << endl;
    // forAll(mesh.C(), cellI)
    // {
    //     porosityField[cellI] = 1.0;
    // }

    const vectorField &centres(mesh.C());   // cell centres
    const scalarField V = mesh.V();         //volume of cells
    vectorField &Usource = UEqn.source();   //source terms

    const vectorField &fluidVelocity(U);

    forAll(structuralElements_memb, Elementi)
    {
        point p0(structuralPositions_memb[structuralElements_memb[Elementi][0]]);
        point p1(structuralPositions_memb[structuralElements_memb[Elementi][1]]);
        point p2(structuralPositions_memb[structuralElements_memb[Elementi][2]]);
        vector e_n(calcNorm(p0,p1,p2));
        
        forAll(centres, cellI)
        {
            if (
                isInPorousZone(centres[cellI], p0, p1, p2) || isInPorous_line(centres[cellI], p0, p1, p2)
               )
            {
                // Info << "DEBUG>>> cellI = "<< cellI << endl;
                // Info << "DEBUG>>> u = "<< fluidVelocity[cellI] << endl;
                scalar theta(calcTheta(p0, p1, p2, fluidVelocity[cellI]));
                // Info << "DEBUG>>> theta = "<< theta << endl;

                scalar cd(calcCd(theta));
                // Info << "DEBUG>>> cd = "<< cd << endl;

                scalar cl(calcCl(theta));
                // Info << "DEBUG>>> cl = "<< cl << endl;

                vector i_d(fluidVelocity[cellI]/mag(fluidVelocity[cellI]));
                // Info << "DEBUG>>> i_d = "<< i_d << endl;

                vector i_l((i_d^e_n)^i_d/mag((i_d^e_n)^i_d));
                // Info << "DEBUG>>> i_l = "<< i_l << endl;

                vector u_oo(sqrt(2/(2-cd-cl))*fluidVelocity[cellI]);
                // Info << "DEBUG>>> u_oo = "<< u_oo << endl;
                
                
                Usource[cellI] -= 0.5*magSqr(u_oo)* V[cellI]/ thickness_memb *(cd*i_d+cl*i_l);
                // Info << ">>> addResistance "<<Usource[cellI] <<"to cell " <<cellI << endl;
                
                Fh[cellI]=0.5*magSqr(u_oo)* V[cellI]/ thickness_memb *(cd*i_d+cl*i_l)*1000.0;
                // porosityField[cellI] = Sn_memb;

            }
        }
    }
}

void Foam::netPanel::creatPoroField(
        volScalarField &porosityField,
        volVectorField &Fh,
        const fvMesh &mesh) const
{
    forAll(mesh.C(), cellI)
        {
            porosityField[cellI] = 1.0;
            Fh[cellI]= vector(0,0,0);
        }
    // get the center of all the cells
    const vectorField &centres(mesh.C());
    //- step 2 assign sn to the proper mesh
    forAll(structuralElements_memb, Elementi)

    {
        point p0(structuralPositions_memb[structuralElements_memb[Elementi][0]]);
        point p1(structuralPositions_memb[structuralElements_memb[Elementi][1]]);
        point p2(structuralPositions_memb[structuralElements_memb[Elementi][2]]);
        forAll(centres, cellI) // loop through all the cell,
        {
            if (isInPorousZone(centres[cellI], p0, p1, p2) || isInPorous_line(centres[cellI], p0, p1, p2))
            {
                porosityField[cellI] = Sn_memb;
                Fh[cellI]=vector(0.5,0,0);
            }
        }
    }
}



// * * * * * * * * * * * * * * Communication Functions  * * * * * * * * * * * * * * //


void Foam::netPanel::readSurf(
        const dictionary &structuralElements)
{
    scalar listLength(readScalar(structuralElements.lookup("numOfSurf")));
    List<vector> surf(listLength, vector::zero);
    forAll(surf, i)
    {
        word surf_name("e" + Foam::name(i));
        structuralElements.lookup(surf_name)>>surf[i] ;
    }
    structuralElements_memb = surf;
    
}


void Foam::netPanel::readPosi(
        const dictionary &structuralPositions)
{
    scalar listLength(readScalar(structuralPositions.lookup("numOfPoint")));
    List<vector> posi(listLength, vector::zero);
    forAll(posi, i)
    {
        word point_name("p" + Foam::name(i));
        structuralPositions.lookup(point_name) >> posi[i];
    }
    structuralPositions_memb = posi;
}



// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

// ************************************************************************* //
