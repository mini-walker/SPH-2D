@echo off
setlocal enabledelayedexpansion

rem ===== User Info: Only set on first-time use =====
rem git config --global user.name "mini-walker"
rem git config --global user.email "sjin@mun.ca"

echo.
echo ===== Step 0: Create .gitignore =====
echo # Auto-generated .gitignore > .gitignore
echo build/>>.gitignore
echo */build/>>.gitignore
echo */amrex_build/>>.gitignore
echo *.o>>.gitignore
echo *.a>>.gitignore
echo *.so>>.gitignore
echo *.exe>>.gitignore
rem echo *.obj>>.gitignore
echo [.gitignore created]

echo.
echo ===== Step 1: Remove all tracked "build" directories =====
for /f "delims=" %%D in ('dir /ad /b /s build 2^>nul') do (
    echo Removing Git tracking from: %%D
    git rm -r --cached "%%D" >nul 2>nul
)
echo [Done]

echo.
echo ===== Step 2: Initialize Git repo if not exists =====
IF NOT EXIST ".git" (
    git init
    git branch -M main
    echo [INFO] Git repository initialized
)

echo.
echo ===== Step 3: Set remote if not already set =====
git remote -v | findstr /C:"origin" >nul
IF ERRORLEVEL 1 (
    git remote add origin git@github.com:mini-walker/SPH-2D.git
    echo [INFO] Remote origin set
) ELSE (
    echo [INFO] Remote origin already exists
)

echo.
echo ===== Step 4: Stage all changes (including deletes) =====
git add -A

echo.
echo ===== Step 5: Auto commit with timestamp =====
for /f %%i in ('powershell -command "Get-Date -Format yyyy-MM-dd_HH:mm:ss"') do set timestamp=%%i
git commit -m "Auto commit at %timestamp%" >nul 2>nul
echo [INFO] Commit done (if there were changes)

echo.
echo ===== Step 6: Pull and merge remote changes if needed =====
git pull origin main --allow-unrelated-histories --no-edit

echo.
echo ===== Step 7: Push to GitHub =====
git push origin main

echo.
echo [SUCCESS] Code has been pushed to GitHub.
pause
