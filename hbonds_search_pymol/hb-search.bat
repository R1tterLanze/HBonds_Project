# Script to start hb-search without setting the environment variable & hb-define file

set "scriptdir=%~dp0"
set PSE_FILE=%scriptdir%\period-table-info.txt
echo %*

"%scriptdir%\Windows\hb-search.exe" -hb "%scriptdir%\hb-define.txt" %*
