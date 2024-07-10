<#
.Description
    This script builds and installs the project dependencies.
.PARAMETER BuildType
    Sets the build type.
.PARAMETER WorkDir
    Path to the directory where packages are downloaded and extracted, repositories are cloned, and targets are built.
.PARAMETER InstallDir
    Path to the directory where the libraries are installed.
.PARAMETER GitTags
    Git tags used to checkout the desired git branches.
    Care must be taken when setting the tag strings. Different repositories use different tag formats.
    Let X, Y and Z be the major, minor and patch release versions, respectively. The tags are formatted as follows:
    zlib       : 'vX.Y.Z'     (https://github.com/madler/zlib/tags)
    hdf5       : 'hdf5-X_Y_Z' (https://github.com/HDFGroup/hdf5/tags)
    sundials   : 'vX.Y.Z';    (https://github.com/LLNL/sundials/tags)
    yamlcpp    : varies       (https://github.com/jbeder/yaml-cpp/tags)
    units      : 'vX.Y.Z'     (https://github.com/LLNL/units/tags)
    catch2     : 'vX.Y.Z'     (https://github.com/catchorg/Catch2/tags)
    googletest : 'vX.Y.Z'     (https://github.com/google/googletest/tags)
    cantera    : 'vX.Y.Z'     (https://github.com/Cantera/cantera/tags)
    fmt        : 'X.Y.Z'      (https://github.com/fmtlib/fmt/tags)
    eigen      : 'X.Y.Z'      (https://gitlab.com/libeigen/eigen/-/tags)
    boost      : 'X.Y.Z'      (https://github.com/boostorg/boost/tags)
.PARAMETER ParallelJobs
    The maximum number of concurrent processes to use when building.
.PARAMETER Clean
    Removes the work directory upon finishing the installation.
.SYNOPSIS
    Used to build the dependencies.
#>

using module ./GenericInstaller.psm1
using module ./BoostInstaller.psm1
using module ./CanteraInstaller.psm1

Param(
    [Parameter(Mandatory = $true)] [string] $WorkDir,
    [Parameter(Mandatory = $true)] [string] $InstallDir,
    [Parameter(Mandatory = $false)] [ValidateSet('Release', 'Debug', 'RelWithDebInfo')] [string] $BuildType = 'Release',
    [Parameter(Mandatory = $false)] [hashtable]$GitTags = @{ `
            zlib       = 'v1.2.13'; `
            hdf5       = 'hdf5-1_14_0'; `
            sundials   = 'v6.7.0'; `
            yamlcpp    = '0.8.0'; `
            units      = 'v0.9.1'; `
            catch2     = 'v3.5.2'; `
            googletest = 'v1.14.0'; `
            cantera    = 'v3.0.0'; `
            fmt        = '10.2.1'; `
            eigen      = '3.4.0'; `
            boost      = '1.84.0'    
    },
    [Parameter(Mandatory = $false)] [ValidateRange(1, [int]::MaxValue)] [int] $ParallelJobs = 4,
    [Parameter(Mandatory = $false)] [Switch] $Clean = $False
)

if (-not(Test-Path -Path $WorkDir)) {
    New-Item $WorkDir -Type Directory
}
$WorkDir = Resolve-Path $WorkDir

$TransriptPath = (Join-Path $WorkDir 'install_log.txt')
Start-Transcript -Path $TransriptPath

Function Invoke-Terminate {
    Stop-Transcript
    Exit
}

# Check version, powershell 7+ is supported
$MajorPSVers = $PSVersionTable.PSVersion.Major
if ([int]$MajorPSVers -lt 7) {
    Write-Error ('Found PowerShell version $MajorPSVers. Version 7+ is required.')
    Invoke-Terminate
    
}

# Check if OS is supported
if ( -not ($IsLinux -or $IsWindows)) {
    Write-Error ('Unsupported operating system. The following are supported: Linux and Windows.')
    Invoke-Terminate
}

if (-not(Test-Path -Path $InstallDir)) {
    New-Item $InstallDir -Type Directory
}
$Global:InstallDir = Resolve-Path $InstallDir

$Global:RepositoryDir = (Join-Path $WorkDir 'repositories')
New-Item -Force $Global:RepositoryDir -Type Directory

$Global:DownloadDir = (Join-Path $WorkDir 'download')
New-Item -Force $Global:DownloadDir -Type Directory

$Global:ExtractDir = (Join-Path $WorkDir 'extract')
New-Item -Force $Global:ExtractDir -Type Directory

$Global:BuildDir = (Join-Path $WorkDir 'build')
New-Item -Force $Global:BuildDir -Type Directory

Write-Host 'Work directory              : ' $WorkDir
Write-Host 'Installation directory      : ' $InstallDir
Write-Host 'Build type                  : ' $BuildType
Write-Host 'Parallel jobs per build     : ' $ParallelJobs
Write-Host 'Tagged of branches to checkout : ' ($GitTags | Out-String)

# =========================================================== zlib

try {
    $zlibInstaller = [GenericInstaller]::new(@{
            Name         = 'zlib'
            Repository   = 'https://github.com/madler/zlib.git'
            Tag          = $GitTags.zlib
            BuildType    = $BuildType
            ParallelJobs = $ParallelJobs
        })

    $zlibInstaller.Run()

    $env:Path += (';' + $zlibInstaller.Directory())

    # =========================================================== hdf5

    $hdf5Installer = [GenericInstaller]::new(@{
            Name         = 'hdf5'
            Repository   = 'https://github.com/HDFGroup/hdf5.git'
            Tag          = $GitTags.hdf5
            BuildType    = $BuildType
            ParallelJobs = $ParallelJobs
            Options      = @(
                '-DBUILD_STATIC_LIBS:BOOL=ON',
                '-DBUILD_SHARED_LIBS:BOOL=OFF',
                '-DHDF5_BUILD_TOOLS:BOOL=OFF',
                '-DHDF5_BUILD_EXAMPLES:BOOL=OFF',
                '-DBUILD_TESTING:BOOL=ON',
                '-DZLIB_USE_EXTERNAL:BOOL=OFF',
                '-DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON',
            ('-DZLIB_ROOT={0}' -f $zlibInstaller.Directory()),
            ('-DZLIB_INCLUDE_DIR:PATH={0}/include' -f $zlibInstaller.Directory()), 
            ('-DZLIB_LIBRARY:FILEPATH={0}/lib/{1}' -f $zlibInstaller.Directory(), ($IsLinux ? "libz.a" : "zlib.lib"))
            )
        })

    $hdf5Installer.Run()

    # =========================================================== sundials

    $sundialsCMakeOptions = @(
        '-DBUILD_STATIC_LIBS:BOOL=ON',
        '-DBUILD_SHARED_LIBS:BOOL=OFF',
        '-DBUILD_CVODE:BOOL=ON',
        '-DBUILD_CVODES:BOOL=ON',
        '-DBUILD_ARKODE:BOOL=OFF',
        '-DBUILD_IDA:BOOL=ON',
        '-DBUILD_IDAS:BOOL=ON',
        '-DBUILD_KINSOL:BOOL=ON',
        '-DBUILD_EXAMPLES:BOOL=ON',
        '-DBUILD_TESTING:BOOL=ON'
    )
    if ($IsLinux) {
        $sundialsCMakeOptions += '-DCMAKE_C_FLAGS=-fPIC'
    }

    $sundialsInstaller = [GenericInstaller]::new(@{
            Name         = 'sundials'
            Repository   = 'https://github.com/LLNL/sundials.git'
            Tag          = $GitTags.sundials
            BuildType    = $BuildType
            ParallelJobs = $ParallelJobs
            Options      = $sundialsCMakeOptions
        })

    $sundialsInstaller.Run()
    
    # =========================================================== yaml

    $yamlcppInstaller = [GenericInstaller]::new(@{
            Name         = 'yaml-cpp'
            Repository   = 'https://github.com/jbeder/yaml-cpp.git'
            Tag          = $GitTags.yamlcpp
            BuildType    = $BuildType
            ParallelJobs = $ParallelJobs
            Options      = @('-DYAML_BUILD_SHARED_LIBS:BOOL=ON')
        })

    $yamlcppInstaller.Run()
    
    # =========================================================== units

    $unitsInstaller = [GenericInstaller]::new(@{
            Name         = 'units'
            Repository   = 'https://github.com/LLNL/UNITS.git'
            Tag          = $GitTags.units
            BuildType    = $BuildType
            ParallelJobs = $ParallelJobs
        })

    $unitsInstaller.Run()
    
    # =========================================================== catch2: included to accommoddate existing tests before migration to google test
    
    $catch2Installer = [GenericInstaller]::new(@{
            Name         = 'catch2'
            Repository   = 'https://github.com/catchorg/Catch2.git'
            Tag          = $GitTags.catch2
            BuildType    = $BuildType
            ParallelJobs = $ParallelJobs
        })

    $catch2Installer.Run()

    # =========================================================== googletest
    $googletestCMakeOptions = @()
    if ($IsWindows) {
        $googletestCMakeOptions = '-Dgtest_force_shared_crt:BOOL=ON'
    }
    
    $googletestInstaller = [GenericInstaller]::new(@{
            Name         = 'googletest'
            Repository   = 'https://github.com/google/googletest.git'
            Tag          = $GitTags.googletest
            BuildType    = $BuildType
            ParallelJobs = $ParallelJobs
            Options      = $googletestCMakeOptions
        })

    $googletestInstaller.Run()

    # =========================================================== fmt
    
    $fmtInstaller = [GenericInstaller]::new(@{
            Name         = 'fmt'
            Repository   = 'https://github.com/fmtlib/fmt.git'
            Tag          = $GitTags.fmt
            BuildType    = $BuildType
            ParallelJobs = $ParallelJobs
            Options      = @(
                '-DCMAKE_POSITION_INDEPENDENT_CODE=TRUE',
                '-DBUILD_SHARED_LIBS=TRUE'
            )
        })

    $fmtInstaller.Run()

    # =========================================================== eigen

    $eigenInstaller = [GenericInstaller]::new(@{
            Name         = 'eigen'
            Repository   = 'https://gitlab.com/libeigen/eigen.git'
            Tag          = $GitTags.eigen
            BuildType    = $BuildType
            ParallelJobs = $ParallelJobs
        })

    $eigenInstaller.Run()

    # =========================================================== boost: uses bootstrap

    $boost = [BoostInstaller]::new(@{
            Tag          = $GitTags.boost
            BuildType    = $BuildType
            ParallelJobs = $ParallelJobs
        })

    $boost.Run()

    # =========================================================== cantera: uses scons

    $DepenendenciesPaths = @{
        "fmt"        = $fmtInstaller.Directory()
        "yamlcpp"    = $yamlcppInstaller.Directory()
        "eigen"      = $eigenInstaller.Directory()
        "googletest" = $googletestInstaller.Directory()
        "sundials"   = $sundialsInstaller.Directory()
        "boost"      = $boost.IncludeDirectory()
    }

    $cantera = [CanteraInstaller]::new(@{
            Tag                 = $GitTags.cantera
            BuildType           = $BuildType
            ParallelJobs        = $ParallelJobs
            DepenendenciesPaths = $DepenendenciesPaths
        })

    $cantera.Run()
}
catch {
    Write-Error ($_.Exception | Format-List -Force | Out-String) -ErrorAction Continue
    Write-Error ($_.InvocationInfo | Format-List -Force | Out-String) -ErrorAction Continue
    Invoke-Terminate
}

Write-Host ( `
        "`nCongratulations! You did it!`n" `
        + "`nYou must now either:`n" `
        + (" - add {0} to your path before building, or`n" -f $Global:InstallDir) `
        + (" - configure the build with -DCMAKE_PREFIX_PATH={0}." -f $Global:InstallDir) `
)

Stop-Transcript

if ($Clean) {
    Remove-Item -Recurse -Force $WorkDir
}