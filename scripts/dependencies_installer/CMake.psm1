class CMake {
    hidden [string] $Name
    hidden [string] $RepositoryDir
    hidden [string] $BuildDir
    hidden [string] $InstallDir
    hidden [string] $BuildType = 'Release'
    hidden [int] $ParallelJobs = 1
    hidden [string[]] $Options = @()
   
    # Default constructor
    CMake() { $this.Initialise(@{}) }

    CMake([hashtable] $Properties) {
        [CMake]::ValidateProperties($Properties)
        $this.Initialise($Properties)
    }

    hidden static [void] ValidateProperties([hashtable] $Properties) {

        [string[]] $requiredProperties = @(
            'Name',
            'RepositoryDir',
            'BuildDir',
            'InstallDir'
        )

        [string[]] $optionalProperties = @(
            'BuildType',
            'ParallelJobs',
            'Options'
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

    hidden [void] Configure() {
        $ConfigureScriptBlock = {
            # param ($RepositoryDir, $BuildDir, $InstallDir, $Options)
            # $output = cmake -S $RepositoryDir -B $BuildDir --install-prefix $InstallDir $Options | Out-String
            # Write-Host $output
            param ($RepositoryDir, $BuildDir, $InstallDir, $Options)
            cmake -S $RepositoryDir -B $BuildDir --install-prefix $InstallDir $Options | Out-Host
        }
        Invoke-Command -ScriptBlock $ConfigureScriptBlock `
            -ArgumentList $this.RepositoryDir, $this.BuildDir, $this.InstallDir, $this.Options
        $ExitCode = $LASTEXITCODE
        if ( $ExitCode -gt 0 ) {
            throw ('Failed to configure the build {0}. Exit code = {1}.' `
                    -f $this.RepositoryDir, $ExitCode.ToString())
        }
    }

    hidden [void] Build() {
        $BuildScriptBlock = {
            param ($BuildDir, $BuildType, $ParallelJobs)
            # $output = cmake --build $BuildDir --config $BuildType --parallel $ParallelJobs | Out-String
            # Write-Host $output
            cmake --build $BuildDir --config $BuildType --parallel $ParallelJobs | Out-Host
        }
        Invoke-Command -ScriptBlock $BuildScriptBlock `
            -ArgumentList $this.BuildDir, $this.BuildType, $this.ParallelJobs
        $ExitCode = $LASTEXITCODE
        if ( $ExitCode -gt 0 ) {
            throw ('Failed to build {0} in {1}. Exit code = {2}.' `
                    -f $this.RepositoryDir, $this.BuildDir, $ExitCode.ToString())
        }
    }

    hidden [void] Install() {
        $InstallScriptBlock = {
            param ($BuildDir)
            # $output = cmake --install $BuildDir | Out-String
            # Write-Host $output
            cmake --install $BuildDir | Out-Host
        }
        Invoke-Command -ScriptBlock $InstallScriptBlock -ArgumentList $this.BuildDir
        $ExitCode = $LASTEXITCODE
        if ( $ExitCode -gt 0 ) {
            throw ('Failed to install {0} in {1}. Exit code = {2}.' `
                    -f $this.BuildDir, $this.InstallDir, $ExitCode.ToString())
        }
    }

    [void] Run() {
        $this.Configure()
        $this.Build()
        $this.Install()
    }
}
