using module ./Git.psm1
using module ./Utilities.psm1

class CanteraInstaller {
    hidden [string] $Tag
    hidden [string] $BuildType
    hidden [int] $ParallelJobs = 1
    hidden [hashtable] $DepenendenciesPaths

    hidden [string] $Name = 'cantera'
    hidden [string] $RepositoryDir = (Join-Path $Global:RepositoryDir $this.Name)
    hidden [string] $BuildDir = (Join-Path $Global:BuildDir $this.Name)
    hidden [string] $InstallDir = (Join-Path $Global:InstallDir $this.Name)
    
    hidden [bool] $DoInstall = $false

    # Default constructor
    CanteraInstaller() { $this.Initialise(@{}) }
    
    CanteraInstaller([hashtable] $Properties) {
        [CanteraInstaller]::ValidateProperties($Properties)
        $this.Initialise($Properties)
        $this.ProcessDirectories()
    }

    hidden static [void] ValidateProperties([hashtable] $Properties) {

        [string[]] $requiredProperties = @(
            'Tag',
            'BuildType'
            'DepenendenciesPaths'
        )

        [string[]] $optionalProperties = @(
            'ParallelJobs'
        )

        $allProperties = $requiredProperties + $optionalProperties

        foreach ($property in $requiredProperties) {
            if (-not $Properties.ContainsKey($property)) {
                throw "Not all required properties are provided"
            }
        }

        foreach ($property in $Properties.Keys) {
            if ($allProperties -notcontains $property) {
                throw "Invalid properties are provided"
            }
        }
    }

    # Initialiser method
    hidden [void] Initialise([hashtable]$Properties) {
        foreach ($Property in $Properties.Keys) {
            $this.$Property = $Properties.$Property
        }
    }

    hidden [void] ProcessDirectories() {
        if (-not(Test-Path -Path $this.InstallDir)) {
            New-Item $this.InstallDir -Type Directory
            if (Test-Path -Path $this.RepositoryDir) {
                Remove-Item $this.RepositoryDir -Recurse -Force
            }
            New-Item $this.RepositoryDir -Type Directory
            if (Test-Path -Path $this.BuildDir) {
                Remove-Item $this.BuildDir -Recurse -Force
            }
            New-Item $this.BuildDir -Type Directory
            $this.DoInstall = $true
        }
    }

    hidden [void] GetSource() {
        $canteraGit = [Git]::new(@{
                Name          = $this.Name
                Repository    = 'https://github.com/Cantera/cantera.git'
                Tag           = $this.Tag
                RepositoryDir = $this.RepositoryDir
            })
        $canteraGit.Run()
    }

    # hidden [void] InstallPythonRequirements() {
    #     python -m pip install scons setuptools wheel ruamel.yaml pytest
    # }

    hidden [void] Build() {
        $ExtraIncDirsArray = @(
            (Join-Path $this.DepenendenciesPaths["fmt"] 'include'),
            (Join-Path $this.DepenendenciesPaths["yamlcpp"] 'include'),
            (Join-Path  $this.DepenendenciesPaths["eigen"] 'include'),
            (Join-Path  $this.DepenendenciesPaths["googletest"] 'include')
        )
        $ExtraIncDirs = ($ExtraIncDirsArray -Join [System.IO.Path]::PathSeparator)
        
        $ExtraLibdirsArray = @(
            (Join-Path  $this.DepenendenciesPaths["fmt"] 'lib'),
            (Join-Path  $this.DepenendenciesPaths["yamlcpp"] 'lib'),
            (Join-Path  $this.DepenendenciesPaths["googletest"] 'lib')
        )
        $ExtraLibDirs = ($ExtraLibDirsArray -Join [System.IO.Path]::PathSeparator)        
        
        $BoostIncDir = $this.DepenendenciesPaths["boost"]
        
        $SundialsInclude = (Join-Path $this.DepenendenciesPaths["sundials"] 'include')
        $SundialsLibdir = (Join-Path $this.DepenendenciesPaths["sundials"] 'lib')

        Push-Location $this.RepositoryDir

        #hdf_support='n' `

        $SConsBuildScriptBlock = {
            scons build `
                --jobs="$($this.ParallelJobs)" ` `
                --directory="$($this.RepositoryDir)" `
                prefix="$($this.InstallDir)" `
                python_package='minimal' `
                system_fmt='y' `
                system_eigen='y' `
                system_yamlcpp='y' `
                system_sundials='y' `
                googletest='system' `
                f90_interface='n' `
                boost_inc_dir="$($BoostIncDir)" `
                sundials_include="$($SundialsInclude)" `
                sundials_libdir="$($SundialsLibdir)" `
                extra_inc_dirs="$($ExtraIncDirs)" `
                extra_lib_dirs="$($ExtraLibDirs)" `
            | Out-Host
        }
        Invoke-Command -ScriptBlock $SConsBuildScriptBlock
        $ExitCode = $LASTEXITCODE
        if ( $ExitCode -ne 0 ) {
            throw ('Failed to build cantera. Exit code = {0}.' -f $ExitCode.ToString())
        }
        Pop-Location
    }

    hidden [void] Install() {
        Push-Location $this.RepositoryDir
        $SConsInstallScriptBlock = {
            scons --directory=$this.RepositoryDir install | Out-Host
        }
        Invoke-Command -ScriptBlock $SConsInstallScriptBlock
        $ExitCode = $LASTEXITCODE
        if ( $ExitCode -ne 0 ) {
            throw ('Failed to install. Exit code = {0}.' -f $ExitCode.ToString())
        }
        Pop-Location
    }

    [string] Directory() { return $this.InstallDir }

    [void] Run() {
        if ($this.DoInstall) {
            Utilities\PrintDecoratedMessage $this.Name "Start"
            $this.GetSource()
            #$this.InstallPythonRequirements()
            $this.Build()
            $this.Install()
            Utilities\PrintDecoratedMessage $this.Name "End"
        }
        else {
            Write-Host $this.Name is already installed.
        }
    }
}
