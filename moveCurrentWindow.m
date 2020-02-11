function moveCurrentWindow()
loadlibrary gdi32.dll Windows.h includepath .
g=calllib('gdi32', 'GetForegroundWindow')
unloadlibrary gdi32
