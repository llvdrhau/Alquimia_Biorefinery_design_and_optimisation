#!/usr/bin/env python
"""
pipeline to reduce, parametrize and balance a biological model in SBML format.

Author: Rik van Rosmalen
"""
from __future__ import division
from __future__ import print_function

import argparse
import logging
import getpass
import os
import os.path

import libsbml
import pint

from databases import Brenda, Equilibrator, Metacyc, Rhea, SabioRK, MetaNetX
from sbmlwrap import Model
from datatracer import Tracer, Task
from identifiers import Identifiers

import balancer


def create_tasks(reaction, organism):
    """Create a task per reaction."""
    data = dict(reaction.identifiers)
    data.update({Identifiers.organism_name: organism,
                 Identifiers.SBMLReactionID: reaction.id,
                 Identifiers.reaction: reaction.name})
    tasks = []

    # For the reactions
    tasks.append(Task(data.copy(), [Identifiers.kinetics_keq, Identifiers.kinetics_vmax,
                                    Identifiers.kinetics_kcat, Identifiers.kinetics_kv]))

    # For the species, reactions combo
    for substrate, count in reaction.reactants:
        d = data.copy()
        d.update(substrate.identifiers)
        d[Identifiers.species] = substrate.name
        d[Identifiers.SBMLSpeciesID] = substrate.id
        tasks.append(Task(d, [Identifiers.kinetics_km, Identifiers.kinetics_kcat,
                              Identifiers.kinetics_kcat_km, Identifiers.kinetics_vmax,
                              Identifiers.kinetics_kv]))
    for modifier in reaction.modifiers:
        # Modifiers. Go for ki or km or nothing?
        # Could be inhibitors, activators, enzymes... ?
        d = data.copy()
        d.update(modifier.identifiers)
        d[Identifiers.species] = modifier.name
        d[Identifiers.SBMLSpeciesID] = modifier.id
        tasks.append(Task(d, [Identifiers.kinetics_ki, Identifiers.kinetics_ka]))
    return tasks


def writeSBtab(filename, tasks, include_mutants=False):
    """Write SBtab parameter output file from finished tasks in Lubitz' format."""
    header = ('QuantityType\tSBMLReactionID\tSBMLSpeciesID\tMean\tStd\t' +
              'Unit\tTemperature\tpH\tMinimum\tMaximum\n')

    # These are all the data types Lubitz' script accept.
    types = {
            Identifiers.kinetics_km: ('Michaelis constant', 'mM'),
            Identifiers.kinetics_ki: ('inhibitory constant', 'mM'),
            Identifiers.kinetics_ka: ('activation constant', 'mM'),
            Identifiers.kinetics_keq: ('equilibrium constant', ''),
            Identifiers.kinetics_kv: ('catalytic rate constant geometric mean', '1/s'),
            Identifiers.kinetics_ce: ('concentration of enzyme', 'mM'),
            Identifiers.kinetics_kcatm: ('product catalytic rate constant', '1/s'),
            Identifiers.kinetics_kcat: ('substrate catalytic rate constant', '1/s'),
            Identifiers.kinetics_vmax: ('forward maximal velocity', 'mM/s'),
            Identifiers.kinetics_vmaxm: ('reverse maximal velocity', 'mM/s'),
            Identifiers.kinetics_a: ('reaction affinity', 'kJ/mol'),
            Identifiers.kinetics_mu0: ('standard chemical potential', 'kJ/mol'),
            Identifiers.kinetics_mu: ('chemical potential', 'kJ/mol'),
            Identifiers.kinetics_c: ('concentration', 'mM')
            }

    # reaction_types = {Identifiers.kinetics_keq, Identifiers.kinetics_kv, Identifiers.kinetics_ce,
    #                   Identifiers.kinetics_kcat, Identifiers.kinetics_kcatm,
    #                   Identifiers.kinetics_vmax, Identifiers.kinetics_vmaxm,
    #                   Identifiers.kinetics_a}
    # species_types = {Identifiers.kinetics_c, Identifiers.kinetics_mu, Identifiers.kinetics_mu0}
    # combined_types = {Identifiers.kinetics_km, Identifiers.kinetics_ki, Identifiers.kinetics_ka}

    ureg = pint.UnitRegistry()
    ureg.define('katal = mol / second')
    ureg.define('M = mol / liter')
    all_count = 0
    all_results = set()
    for task in tasks:
        SBMLReactionID = task.data.get(Identifiers.SBMLReactionID, None)
        SBMLSpeciesID = task.data.get(Identifiers.SBMLSpeciesID, None)
        for result in task.kinetic_data:
            all_count += 1
            try:
                type_id, unit = types[result.kinetics_type]
            except KeyError:
                logging.info('Result ignored because of incompatible kinetic constant type.')
                # print('Result ignored because of incompatible kinetic constant type.')
                # print(result.source_query)
                # print(result.source_pubmed or result.source_publication)
                # print(result._sabio)
                continue
            if unit == '' and result.kinetics_unit == '-':
                mean = float(result.kinetics_value)
            elif result.kinetics_unit is None or not result.kinetics_unit:
                logging.info('Result ignored because of missing units.')
                # print('Result ignored because of missing units.')
                # print(result.source_query)
                # print(result.source_pubmed or result.source_publication)
                # print(result._sabio)
                continue
            elif result.mutant and not include_mutants:
                logging.info('Result ignored because of mutant phenotype.')
                # print('Result ignored because of mutant phenotype.')
                # print(result.source_query)
                # print(result.source_pubmed or result.source_publication)
                # print(result._sabio)
                continue
            else:
                try:
                    mean = convertUnits(float(result.kinetics_value), result.kinetics_unit,
                                        unit, ureg)
                except pint.DimensionalityError:
                    logging.info('Result ignored because of unit mismatch.')
                    # print('Result ignored because of unit mismatch.')
                    # print('\t', repr(result.kinetics_value))
                    # print('\t', repr(result.kinetics_unit))
                    # print('\t', repr(unit))
                    # print(result.source_query)
                    # print(result.source_pubmed or result.source_publication)
                    # print(result._sabio)
                    continue
                except Exception:
                    print(repr(result.kinetics_value))
                    print(repr(result.kinetics_unit))
                    print(repr(unit))
                    raise
            std = result.kinetics_std
            if std:
                std = convertUnits(float(result.kinetics_std), result.kinetics_unit,
                                   unit, ureg)
            source = result.source_pubmed or result.source_publication
            s = '\t'.join(str(i) for i in (type_id,
                                           SBMLReactionID,
                                           SBMLSpeciesID,
                                           mean,
                                           std,
                                           result.kinetics_temp,
                                           result.kinetics_ph,
                                           None,
                                           None,
                                           ' '.join((result.source_name, str(source))))) + '\n'
            all_results.add((s, source))

    logging.info("Used {} out of {} results.".format(len(all_results), all_count))
    with open(filename, 'w') as outfile:
        outfile.write(header)
        for result, source in sorted(all_results):
            outfile.write(result)


def convertUnits(ammount, unit_from, unit_to, ureg=None):
    """Convert amount from unit to another unit."""
    if ureg is None:
        ureg = pint.UnitRegistry()
        ureg.define('katal = mol / second')
        ureg.define('M = mol / liter')
    return (ammount * ureg(unit_from)).to(unit_to).magnitude


def write_json(filename, tasks):
    """Write the tasks results to a JSON file."""
    # Recursive encoder snippet from tobique @ http://stackoverflow.com/a/35483750.
    import json
    import inspect

    exclude = {'path', 'wanted', 'available', 'report'}

    class ObjectEncoder(json.JSONEncoder):
        def default(self, obj):
            if hasattr(obj, "to_json"):
                return self.default(obj.to_json())
            elif hasattr(obj, "__dict__"):
                d = dict(
                    (key, value)
                    for key, value in inspect.getmembers(obj)
                    if not key.startswith("__")
                    and not inspect.isabstract(value)
                    and not inspect.isbuiltin(value)
                    and not inspect.isfunction(value)
                    and not inspect.isgenerator(value)
                    and not inspect.isgeneratorfunction(value)
                    and not inspect.ismethod(value)
                    and not inspect.ismethoddescriptor(value)
                    and not inspect.isroutine(value)
                    and key not in exclude
                )
                return self.default(d)
            return obj

    # TODO: Fix possible unicoderror (test with iTO977 model)
    try:
        with open(filename, 'w') as outfile:
            return json.dump(tasks, outfile, cls=ObjectEncoder, indent=2, sort_keys=True)
    except UnicodeEncodeError:
        print("UnicodeEncodeError is JSON dump. Dump aborted")


def checkModifiers(model, default_sbo=None):
    """Check all the reaction modifiers and assign an SBO term.

    If provided, default_sbo will be used. Else, the user will be prompted."""
    for reaction in model.reactions.values():
        for modifier in reaction.modifiers:
            # Check SBO terms.
            sbo = modifier.sbml_object.getSBOTerm()
            inhibition_sbos = [20, 206, 207, 536, 537]  # These are the SBO term Lubitz
            activation_sbos = [13,  21, 459, 461, 462]  # accepts for balancing
            # -1 means SBO term was not set.
            if sbo == -1 or (sbo not in inhibition_sbos and sbo not in activation_sbos):
                # Activation or inhibition?
                if default_sbo:
                    sbo = default_sbo
                else:
                    # No default. Prompt user?
                    while True:
                        try:
                            print("Please choose an SBO term (integer)")
                            print(reaction)
                            print(modifier)
                            sbo = int(raw_input('SBO term: '))
                            break
                        except Exception as e:
                            print(e)
                            print("ERROR: Please try again!")
                # Update models.
                if sbo == -1:
                    modifier.sbml_object.setSBOTerm(sbo)


def runBalancing(modelfile, SBtabfile, outputfile):
    """Run the parameter balancing algorithm from Timo Lubtiz.

    Arguments:
    modelfile -- Location of the SBML model file.
    SBtabfile -- Location of the SBtab file with prior data.
    outputfile -- Location to write the SBtab file with balanced data.
    """
    tl.startTheBusiness(modelfile, SBtabfile, outputfile)


def runKineticizer(modelfile, SBtabfile, outputfile, args=None):
    """Run the kineticizer to add rate laws and parameters to the model. (Lubitz' method).

    Arguments:
    modelfile -- Location of the SBML model file.
    SBtabfile -- Location of the SBtab file with prior data.
    outputfile -- Location to write the SBtab file with balanced data.
    """
    sbml_model = libsbml.readSBML(modelfile)

    with open(SBtabfile, 'r') as inputfile:
        sbtab = inputfile.read().split('\n')

    sbtab = SBtab.SBtabTable(sbtab, SBtabfile, 'QuantityType')

    kineticizer.kineticizer_cs(sbml_model.getModel(), sbtab, 'hal', False,
                               'complete_inh', 'complete_act', True)

    with open(outputfile, 'w') as outfile:
        outfile.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        outfile.write(sbml_model.toSBML())


def getTracerResultStats(model, tracer):
    """Gather the results on how well the data gathering performed."""
    raise NotImplementedError
    possible = {'kinetics_Keq', 'kinetics_Vmax+', 'kinetics_Km', 'kinetics_kcat+',
                'kinetics_kcat_km', 'kinetics_Ki', 'Kinetics_Ka'}
    # Dictionary of (parameter, reaction, substrate): measurements
    results = {}
    # Possible value to find.
    for reaction in model.reactions.values():
        results[('kinetics_Keq', reaction.id, None)] = []
        results[('kinetics_Vmax+', reaction.id, None)] = []
        results[('kinetics_kcat+', reaction.id, None)] = []
        for substrate, count in reaction.reactants + reaction.products:
            results[('kinetics_Km', reaction.id, substrate.id)] = []
            results[('kinetics_kcat+', reaction.id, substrate.id)] = []
            results[('kinetics_kcat_km', reaction.id, substrate.id)] = []
        for modifier in reaction.modifiers:
            results[('kinetics_Ki', reaction.id, modifier.id)] = []
            results[('kinetics_Ka', reaction.id, modifier.id)] = []
    # Special case for concentration since this comes from the sbml directly if found.
    for species in model.species.values():
        results[('concentration', None, species.id)] = []
        if species.initial_concentration >= 0:
            results[('concentration', None, species.id)].append(species.initial_concentration)

    # Obtained
    for task in tracer.tasks:
        reaction, species = None, None
        if Identifiers.SBMLSpeciesID in task.data:
            species = task.data[Identifiers.SBMLSpeciesID]
        if Identifiers.SBMLReactionID in task.data:
            reaction = task.data[Identifiers.SBMLReactionID]
        for parameter in task.available & possible:
            if (parameter, reaction, species) not in results:
                print("{} found but not used".format((parameter, reaction, species)))
            else:
                results[(parameter, reaction, species)].append(task.data[parameter])

    return results


def main(modelpath, organism, outdir, suffix='', databases='all', method='default',
         stopat=False, skipbrendalogin=False, task_dump=False):
    """Run the pipeline."""
    databases = databases.lower()
    logging.info("Loading model")
    # Set up model and create tasks
    m = Model(modelpath)

    logging.info("Creating Tasks")

    tasks = []
    for reaction in m.reactions.values():
        t = create_tasks(reaction, organism)
        if t:
            tasks.extend(t)
            # if len(tasks) > 10:
            #     break

    # Set up our data sources
    transforms = []
    if 'brenda' in databases or databases == 'all':
        logging.info("Brenda: Username/Password required to access Brenda database.")
        if skipbrendalogin:
            # Usefull if we've already cached the whole run previously.
            brenda = Brenda('email', 'password')
        else:
            brenda = Brenda(raw_input("Enter Brenda username: "),
                            getpass.getpass("Enter Brenda password: "))
        transforms.extend(brenda.getTransforms())
    if 'equilibrator' in databases or databases == 'all':
        equilibrator = Equilibrator('data/local/kegg_reactions_CC_ph7.0.csv')
        transforms.extend(equilibrator.getTransforms())
    if 'rhea' in databases or databases == 'all':
        rhea = Rhea('data/local/rhea2xrefs.tsv')
        transforms.extend(rhea.getTransforms())
    if 'metacyc' in databases or databases == 'all':
        metacyc = Metacyc()
        transforms.extend(metacyc.getTransforms())
    if 'sabiork' in databases or databases == 'all':
        sabio = SabioRK()
        transforms.extend(sabio.getTransforms())
    if 'metanetx' in databases or databases == 'all':
        metanetx = MetaNetX()
        transforms.extend(metanetx.getTransforms())

    tracer = Tracer(tasks, transforms)

    # Gather our data
    tracer.run()

    # Dump json task trace if wanted.
    if task_dump:
        write_json(os.path.join(outdir, 'tasks.json'), tracer.tasks)

    # Save the data gathered
    parameterfile = os.path.join(outdir, 'parameters{}.tsv'.format(suffix))
    logging.info("TRACING: Writing found parameters to file ({})".format(parameterfile))
    writeSBtab(parameterfile, tasks)
    results = None  # getTracerResultStats(m, tracer)

    if method == "Lubitz":
        logging.info("PIPELINE: Using Lubitz' balancing and kineticizing.")

        # Import Lubitz' method. This is not installed by default.
        try:
            from lubitz_balancer import SBtab
            from lubitz_balancer import kineticizer
            from lubitz_balancer import balance_reaction_parameters as tl
        except ImportError:
            logging.error("Lubitz' method could not be used because it is not installed.")
            raise

        if stopat == 'check_sbo':
            return results
        # Check modifiers
        logging.info("MODEL: Checking for unassigned reaction modifiers.")
        checkModifiers(m, -1)

        sbomodelfile = os.path.join(outdir, 'model_with_sbos{}.xml'.format(suffix))
        logging.info("MODEL: Added SBO terms to modifiers. Saved in {}".format(sbomodelfile))
        m.writeSBML(sbomodelfile)

        if stopat == 'balance':
            return results
        # Run the parameter balancing
        logging.info("BALANCING: Started parameter balancing.")

        balancedfile = os.path.join(outdir, 'parameters_balanced_lubitz{}.tsv'.format(suffix))
        runBalancing(sbomodelfile, parameterfile, balancedfile)

        logging.info("BALANCING: Finished parameter balancing.")

        if stopat == 'kineticize':
            return results
        # # Create the new model (Warning: Slow!)
        logging.info("MODEL: Started creating kinetic model.")

        kineticmodelfile = os.path.join(outdir, 'kinetic_model_lubitz{}.xml'.format(suffix))
        runKineticizer(sbomodelfile, balancedfile, kineticmodelfile)

        logging.info("MODEL: Finished creating kinetic model.")

        return results
    else:
        logging.info("PIPELINE: Using default balancing and kineticizing.")

        if stopat == 'balance':
            return results
        # Run the parameter balancing
        logging.info("BALANCING: Started parameter balancing.")

        balancedfile = os.path.join(outdir, 'parameters_balanced{}.tsv'.format(suffix))
        # runBalancing(sbomodelfile, parameterfile, balancedfile)

        # This should be wrapped up nicer....
        parameters = {'equilibrium constant': 'keq',
                      'substrate catalytic rate constant': 'kcat',
                      'Michaelis constant': 'km',
                      'concentration': 'c',
                      'inhibitory constant': 'ki',
                      'forward maximal velocity': 'vmax',
                      'activation constant': 'ka',
                      }

        S = m.get_stoichiometry_matrix()
        r_pos = dict(((j, i) for i, j in enumerate(m.reactions.keys())))
        s_pos = dict(((j, i)for i, j in enumerate(m.species.keys())))

        compounds = list(m.species.keys())
        reactions = list(m.reactions.keys())

        data = []
        with open(parameterfile, 'r') as datafile:
            datafile.readline()  # Skip header
            for line in datafile:
                line = line.split('\t')
                if line[3] == 'nan':
                    continue
                v = float(line[3])
                if v < 1E-12:
                    v = 1E-12
                sd = 0.1 * v if line[4] == 'None' else float(line[4])
                p = parameters[line[0]]
                c = None if (line[2] == 'None' or line[2] == 'nan') else line[2]
                r = None if (line[1] == 'None' or line[1] == 'nan') else line[1]
                if (c in s_pos or c is None) and (r in r_pos or r is None):
                    data.append((v, sd, p, c, r))

        result = balancer.balance(compounds, reactions, S, data,
                                  balancer.priors, balancer.dependencies,
                                  balancer.nonlog, balancer.R, balancer.T,
                                  True, s_pos, r_pos)
        x_post = result.median
        sigma_post = result.sd
        x_post_columns = result.columns
        balancer.save_balancing_results(x_post, sigma_post, x_post_columns, balancedfile)
        ###
        logging.info("BALANCING: Finished parameter balancing.")

        if stopat == 'kineticize':
            return results
        # # Create the new model (Warning: Slow!)
        logging.info("MODEL: Started creating kinetic model.")

        kineticmodelfile = os.path.join(outdir, 'kinetic_model{}.xml'.format(suffix))
        m.kineticize(balancedfile)
        m.writeSBML(kineticmodelfile)

        logging.info("MODEL: Finished creating kinetic model.")

        return results


if __name__ == "__main__":
    # Set up command line arguments
    parser = argparse.ArgumentParser(description='Pipeline for model parameterization.')
    parser.add_argument('modelpath',
                        help="path to model file (as sbml)")
    parser.add_argument('organism',
                        help='Organism name. Make sure to match Brenda names.')
    parser.add_argument('outdir', help="path to dir to be created for output files.")
    parser.add_argument('-lubitz', help="Use Lubitz' method for balancing and creating the"
                        " kinetic model", action='store_true')
    parser.add_argument('-v', help='verbose', action='store_true')
    parser.add_argument('-log', help="log to file")

    args = parser.parse_args()

    # Set up logging
    if args.v and args.log:
        logging.basicConfig(level=logging.INFO, filename=args.log, filemode='w')
    elif args.v:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING)

    if args.lubitz:
        method = 'Lubitz'
    else:
        method = 'default'

    # Set up output directory
    try:
        os.makedirs(args.outdir)
    except OSError:
        if not os.path.isdir(args.outdir):
            raise

    # Standard only brenda and kegg are used. The others require some setup [See README]
    main(args.modelpath, args.organism, args.outdir, suffix='',
         databases='brenda kegg', method=method, stopat=False, skipbrendalogin=False)


# Example Hynne glycolis model
# python pipeline.py models/Hynne2001.xml "Saccharomyces cerevisiae" output-hynne/ -v

# Example Costa Lactococcus model
# python pipeline.py models/Costa2014.xml "Lactococcus lactis" output-costa/ -v
