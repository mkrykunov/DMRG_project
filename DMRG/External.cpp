#include <iostream>
#include <process.h>
#include <windows.h>


char MessageString[1024];

void SetRange(int arg1, int arg2)
{
}

int SetPos(int arg)
{
    return arg;
}

int SetStep(int arg)
{
    return arg;
}

void StepIt()
{
}

void Message(const char* msg)
{
    MessageBeep(0xFFFFFFFF);
    std::cout << msg << std::endl;
}

void Error(const char* msg)
{
    MessageBeep(0xFFFFFFFF);
    std::cerr << msg << std::endl;
    exit(1);
}

void SetTrace()
{
}

void Trace(const char* msg)
{
    std::cout << msg << std::endl;
}

void RunProgram(const char* CmdLine)
{
    WinExec(CmdLine, SW_SHOWNORMAL);
}
