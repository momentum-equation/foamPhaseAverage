/*---------------------------------------------------------------------------*\
 *=========                 |
 *\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
 *\\     /   O peration     |
 *\\    /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
 *\\   /     M anipulation  |
 * -------------------------------------------------------------------------------
 license
 This file is part of OpenFOAM.

 OpenFOAM is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2 of the License, or (at your
 option) any later version.

 OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 for more details.

 You should have received a copy of the GNU General Public License
 along with OpenFOAM; if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
   foamPhaseAverage

Author
  Suraj Kulkarni
  email: surajkoolkarni@yahoo.com

Description
Calculates phase average of the specified volField

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "in_out_helpers.H"

template<class FieldType>
void calcPhaseAverage
(
    fvMesh& mesh,
    IOobject& fieldHeader,
    const word& fieldName,
    Time& runTime,
    instantList& timeDirs,
    bool& done,
    scalar& phaseStartTime,
    scalar& cycleTime
);

int main(int argc, char *argv[])
{
    Foam::timeSelector::addOptions();

    #include "addRegionOption.H"

    Foam::argList::validArgs.append("fieldName");

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    runTime.setTime(timeDirs[0], 0);

    #include "createNamedMesh.H"

    // get filename from command line
    word fieldName = args[1];

    word dictionaryFileName = "phaseAverageDict";
    word dictionaryDirName = "constant";

    IOobject phaseAverageDictIO = generateIOobject(dictionaryFileName,mesh,dictionaryDirName);

    // Check the if the dictionary is present and follows the OF format
    if (!phaseAverageDictIO.headerOk())
        FatalErrorIn(args.executable()) << "Cannot open specified dictionary "
            << dictionaryFileName << exit(FatalError);


    const dictionary phaseAverageDict = IOdictionary(phaseAverageDictIO);

    scalar phaseStartTime(readScalar(phaseAverageDict.lookup("phaseStartTime")));
    scalar cycleTime(readScalar(phaseAverageDict.lookup("cycleTime")));

    bool done = false;

    IOobject fieldHeader
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (fieldHeader.headerOk())
    {
        calcPhaseAverage<volScalarField>(mesh, fieldHeader, fieldName, runTime, timeDirs, done, phaseStartTime, cycleTime);
        calcPhaseAverage<volVectorField>(mesh, fieldHeader, fieldName, runTime, timeDirs, done, phaseStartTime, cycleTime);
        calcPhaseAverage<volTensorField>(mesh, fieldHeader, fieldName, runTime, timeDirs, done, phaseStartTime, cycleTime);
        calcPhaseAverage<volSymmTensorField>(mesh, fieldHeader, fieldName, runTime, timeDirs, done, phaseStartTime, cycleTime);
        calcPhaseAverage<volSphericalTensorField>(mesh, fieldHeader, fieldName, runTime, timeDirs, done, phaseStartTime, cycleTime);
    }

    Info<< "End\n" << endl;

    return 0;
}

template<class FieldType>
void calcPhaseAverage
(
    fvMesh& mesh,
    IOobject& fieldHeader,
    const word& fieldName,
    Time& runTime,
    instantList& timeDirs,
    bool& done,
    scalar& phaseStartTime,
    scalar& cycleTime
)
{
    label nfield = 0;
    const word meanFieldName = fieldName + "_phase_locked_" + name(phaseStartTime);

    if(fieldHeader.headerOk())
    {
      if(!done && fieldHeader.headerClassName()== FieldType::typeName)
      {
          Info << "meanField = " << meanFieldName << endl;

          FieldType dummy
          (
              IOobject
              (
                  fieldName,
                  runTime.timeName(),
                  mesh,
                  IOobject::MUST_READ
              ),
              mesh
          );

          FieldType meanField
          (
              IOobject
              (
                  meanFieldName,
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ
              ),
              dummy
          );

          meanField *= scalar(0.0);

          for (label timeI = 0; timeI < timeDirs.size(); timeI++)
          {
            runTime.setTime(timeDirs[timeI], timeI);

            if(timeDirs[timeI] == Foam::instant(phaseStartTime))
            {
                Info << "Time = " << runTime.timeName() <<endl;

                IOobject io
                (
                    fieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ
                );

                if (io.headerOk())
                {
                    mesh.readUpdate();

                    if(!done && io.headerClassName() == FieldType::typeName)
                    {
                        Info << "   Reading " << io.headerClassName() << " " <<io.name() << endl;

                        FieldType field(io, mesh);

                        meanField += field;
                        nfield++;
                    }
                }
                else
                {
                    Info << "   No Field " << fieldName << endl;
                }

                phaseStartTime += cycleTime;
            }
          }

          if(nfield > 0)
          {
              Info << "number of fields added = " << nfield << endl;
              meanField /= nfield;
          }

          Info<< "writing to timeDir " << runTime.timeName()  << endl;
          meanField.write();
          done = true;
        }
    }

}




// ************************************************************************* //
