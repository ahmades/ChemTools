class Git {
    hidden [string] $Name
    hidden [string] $Repository
    hidden [string] $Tag
    hidden [string] $RepositoryDir
      
    # Default constructor
    Git() { $this.Initialise(@{}) }

    Git([hashtable] $Properties) {
        [Git]::ValidateProperties($Properties)
        $this.Initialise($Properties)
    }

    # Initialiser method
    hidden [void] Initialise([hashtable]$Properties) {
        foreach ($Property in $Properties.Keys) {
            $this.$Property = $Properties.$Property
        }
    }

    hidden static [void] ValidateProperties([hashtable] $Properties) {

        [string[]] $requiredProperties = @(
            'Name', 
            'Repository', 
            'Tag'
            'RepositoryDir'
        )

        foreach ($property in $requiredProperties) {
            if (-not $Properties.ContainsKey($property)) {
                throw "Not all required properties are provided"
            }
        }
    }
    
    hidden [bool] Clone() {
        # Clone fails if the repo has been already cloned. Try to pull first then fallback to cloning if pulling fails.
        $PullScriptBlock = {
            param ($RepositoryDir, $Tag)
            git -C $RepositoryDir pull origin refs/tags/$Tag | Out-Host
        }
        Invoke-Command -ScriptBlock $PullScriptBlock -ArgumentList $this.RepositoryDir, $this.Tag
        if ( $LASTEXITCODE -eq 0 ) {
            # repo exists and the correct tag is checked out, nothing to do
            return $false
        }
        
        $CloneScriptBlock = {
            Param ($Repository, $RepositoryDir)
            git clone $Repository $RepositoryDir | Out-Host
            git -C $RepositoryDir status | Out-Host
        }
        Invoke-Command -ScriptBlock $CloneScriptBlock -ArgumentList $this.Repository, $this.RepositoryDir
        
        $ExitCode = $LASTEXITCODE
        if ( $LASTEXITCODE -gt 0 ) {
            throw ('Failed to clone {0} to {1}. Exit code = {2}.' `
                    -f $this.Repository, $this.RepositoryDir, $ExitCode.ToString())
        }
        return $true
    }

    hidden [void] Checkout() {
        # fetch all tags
        $FetchTagsScriptBlock = {
            param ($RepositoryDir)
            git -C $RepositoryDir fetch --all --tags | Out-Host
            git -C $RepositoryDir status | Out-Host
        }
        Invoke-Command -ScriptBlock $FetchTagsScriptBlock -ArgumentList $this.RepositoryDir
        $ExitCode = $LASTEXITCODE
        if ( $ExitCode -gt 0 ) {
            throw ('Failed to fetch all tags from {0}. Exit code = {1}.' `
                    -f $this.Repository, $ExitCode.ToString())
        }
    
        # Check out the branch with the required tag
        $CheckOutTaggedBranchScriptBlock = {
            param($RepositoryDir, $Tag)
            git -C $RepositoryDir checkout tags/$Tag -b $Tag | Out-Host
            git -C $RepositoryDir status | Out-Host
        }
        Invoke-Command -ScriptBlock $CheckOutTaggedBranchScriptBlock -ArgumentList $this.RepositoryDir, $this.Tag
        $ExitCode = $LASTEXITCODE
        if ( $ExitCode -gt 0 ) {
            throw ('Failed to checkout tags/{0} from {1}. Exit code = {2}.' `
                    -f $this.Tag, $this.Repository, $ExitCode.ToString())
        }
    }

    [void] Run() {       
        if ($this.Clone()) {
            $this.Checkout()
        }
    }
}
