Info<< "Reading mass flow properties\n" << endl;

IOdictionary massFlowDict
(
    IOobject
    (
        "massFlowProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

wordList inletPatches(massFlowDict.lookup("inletPatches"));
wordList outletPatches(massFlowDict.lookup("outletPatches"));

fileName path = runTime.caseName();

OFstream outInfo ("outInfo.case");

outInfo << "time   "
        << "m_in   "
        << "v_in   "
        << "H_in   "
        << "m_out  "
        << "v_out  "
        << "H_out  "
        << "m_net  "
        << "v_net  "
        << "H_net  "
        << endl;
