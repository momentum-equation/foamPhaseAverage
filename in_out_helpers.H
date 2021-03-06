template<class FieldType>
bool headerChecker
(
    IOobject& fieldHeader
)
{
    bool check = false;
    check = fieldHeader.typeHeaderOk<FieldType>(true);
    return check;
}

Foam::IOobject generateIOobject
( 
    word& dictionaryName,
    fvMesh& mesh,
    word& dictionaryLoc   
)
{

    fileName dictionaryDir = mesh.time().path()/dictionaryLoc;
    
    if(!isDir(dictionaryDir.path()))
    {
        mkDir(dictionaryDir);
    }

    IOobject dictIO
    (
        dictionaryName,
        dictionaryDir,
        mesh,
        IOobject::MUST_READ
    );
    
    return dictIO;
}  
