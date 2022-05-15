# ---------------------------- #
#  Important Shared Functions  #
# ---------------------------- #

#def get_active_cohorts():
#    """Returns a list of active cohorts.
#    Active cohorts are defined in config.yml data > cohorts > [cohortname] > active: True
#    """
#    return(c for c in config['data']['cohorts'] if config['data']['cohorts'][c]['active'])

def get_cohort_config(cohort):
    """Returns the configuration details for specific cohort.
    Importantly, retrieving a cohort using this function (rather than directly
    accessing it via the config dictionary) will also inject the default settings
    (config > data > defaults) wherever a specific setting is not specified.
    """
    config_data = dict(config['data']['defaults'])
    if 'settings' in config['data']['cohorts'][cohort]:
        for key in config['data']['cohorts'][cohort]['settings']:
            if key in config_data:
                config_data[key] = config['data']['cohorts'][cohort]['settings'][key]
            else:
                raise(Exception('Setting {} does not exist for cohort {}. Available settings: {}'.format(
                    key, cohort, ', '.join(config_data.keys())
                )))
    return(config_data)

def get_cohort_data(cohort):
    """Parses the samplesheet for a specific cohort.
    """
    samplesheet = pd.read_csv(config['data']['cohorts'][cohort]['samplesheet'], comment='#').drop_duplicates()
    return(samplesheet)

def get_all_samples(cohort=None):
    """Retrieves all samples to be processed.
    Does so by calling get_cohort_data

    Keyword arguments:
        cohort -- Name of a cohort, OPTIONAL. If not specified, returns all samples
                  across all cohorts.
    """
    all_samples = pd.concat([
        get_cohort_data(cohort_name).assign(cohort_name = cohort_name)
        for cohort_name
        in config['data']['cohorts']
        if config['data']['cohorts'][cohort_name]['active']
        ])

    if cohort is None:
        return(all_samples)
    else:
        return(all_samples[all_samples.cohort_name == cohort])

#def get_all_samples_list():
#    """Returns all samples in list format.
#    By calling get_all_samples()
#    """
#    return(get_all_samples().sample_name.unique().tolist())

def get_all_samples_with_cohorts():
    """Returns a list of tuples with cohort name and sample name."""
    samples = get_all_samples()[['cohort_name', 'sample_name']]
    return(zip(
        samples.drop_duplicates().cohort_name.tolist(),
        samples.drop_duplicates().sample_name.tolist()
    ))

def clean(command):
    """Cleans a snakemake command string by replacing whitespace with a single space.
    Useful for stripping down multiline shell commands so that they read more
    cleanly on snakemake's printed output.
    """
    command = command.replace('\n', ' ').replace('\t', ' ')
    while '  ' in command:
        command = command.replace('  ', ' ')
    return(command.strip())

