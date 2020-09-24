/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

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
    Foam::argList::validArgs.append("fieldType");

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    runTime.setTime(timeDirs[0], 0);

    #include "createNamedMesh.H"

    // get filename from command line
    word fieldName = args[1];
    word fieldType = args[2];

    word dictionaryFileName = "phaseAverageDict";
    word dictionaryDirName = "constant";

    IOobject phaseAverageDictIO = generateIOobject(dictionaryFileName,mesh,dictionaryDirName);

    // Check the if the dictionary is present and follows the OF format
    if (!phaseAverageDictIO.typeHeaderOk<dictionary>(true))
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

    if (fieldType=="scalar")
    {
        calcPhaseAverage<volScalarField>(mesh, fieldHeader, fieldName, runTime, timeDirs, done, phaseStartTime, cycleTime);
    }
    else if (fieldType =="vector")
    {
        calcPhaseAverage<volVectorField>(mesh, fieldHeader, fieldName, runTime, timeDirs, done, phaseStartTime, cycleTime);
    }
    else if (fieldType == "tensor")
    {
        calcPhaseAverage<volTensorField>(mesh, fieldHeader, fieldName, runTime, timeDirs, done, phaseStartTime, cycleTime);
    }
    else if (fieldType == "symmTensor")
    {
        calcPhaseAverage<volSymmTensorField>(mesh, fieldHeader, fieldName, runTime, timeDirs, done, phaseStartTime, cycleTime);
    }
    else if (fieldType == "sphTensor")
    {
        calcPhaseAverage<volSphericalTensorField>(mesh, fieldHeader, fieldName, runTime, timeDirs, done, phaseStartTime, cycleTime);
    }

    else
    {
        FatalErrorIn(args.executable()) << "Cannot find specified fieldType: "
            << fieldType << exit(FatalError);
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

    if(headerChecker<FieldType>(fieldHeader))
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

                if (headerChecker<FieldType>(io))
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
    else
    {

    }

}




// ************************************************************************* //
