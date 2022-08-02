/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
 No.
Class
    Foam::netPanel

Description

SourceFiles
    netPanel.H
    netPanel.C

\*---------------------------------------------------------------------------*/
#ifndef netPanel_H
#define netPanel_H
// class guards that prevents double declaration

#include "dictionary.H"
#include "cellSet.H"
#include "fvMesh.H"
#include "scalar.H"
#include "fvMatrices.H"

// might be used
// #include "pointField.H"
// #include "volFields.H"
// #include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class netPanel Declaration
\*---------------------------------------------------------------------------*/

class netPanel
{
protected:
    //- Private data
    const dictionary netDict_memb;
    const scalar Sn_memb;        // solidity of the netting
    const scalar thickness_memb; // thickness of the net panel in CFD, thickness of the porous media.
    const scalar dw_memb;        // diameter of the twine in netting
    
    const scalar ropeEnhance_memb;    // 0: no enhance, might be a hole;
                                      // 1: the thickness of twine is the same with thickness of net;
                                      // 2: enhanced, the thickness of twine is twice the thickness of net.


    List<vector> structuralElements_memb;        // screen model of netting
    List<vector> structuralPositions_memb;       // positions of each noeds
    

    //-Private member function
    vector calcNorm(
        const point &pointI,
        const point &pointII,
        const point &pointIII) const;

    scalar calcTheta(
            const point& pointI,
            const point& pointII,
            const point& pointIII,
            const vector& fluidVelocity)const;
    
    scalar calcArea(
        const point &pointI,
        const point &pointII,
        const point &pointIII) const;
    
    scalar calcDist(
        const point &pointI,
        const point &pointII) const;
    
    scalar calcDistanceFromPoint2Panel(
            const point &x,
            const vector &structuralElementi) const;

    bool isInPorous_line(
            const point &x,
            const vector &point_a,
            const vector &point_b,
            const vector &point_c) const;

    bool isInPorousZone(
            const point &x,
            const vector &point_a,
            const vector &point_b,
            const vector &point_c) const;


    scalar calcCd(const scalar &theta) const
    {
        return 0.04+cos(theta)*(-0.04+Sn_memb-1.24*pow(Sn_memb,2)+13.7*pow(Sn_memb,3));
    };
    scalar calcCl(const scalar &theta) const
    {
        return sin(2*theta)*(0.57*Sn_memb-3.54 * pow(Sn_memb, 2) + 10.1 * pow(Sn_memb, 3));
    };


public:
    // Constructors
    netPanel( const dictionary &netDict);

    //- Destructor
    ~netPanel();

    // public Member Functions
    void addResistance(
        const volVectorField &U,
        volVectorField &Fh,
        fvVectorMatrix &UEqn,
        // volScalarField &porosityField,
        const fvMesh &mesh) const;

    void creatPoroField(
        volScalarField &porosityField,
        volVectorField &Fh,
        const fvMesh &mesh) const;


    scalar porosity() const
    {
        return 1-Sn_memb*dw_memb/thickness_memb*atan(1.0);
    };

    scalar Sn() const
    {
        return Sn_memb;
    };



// * * * * * * * * * * * * * * Communication Functions  * * * * * * * * * * * * * * //

    void readPosi(const dictionary &structuralPositions);
    void readSurf(const dictionary &structuralElements);
    
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "netPanelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //r&  t
