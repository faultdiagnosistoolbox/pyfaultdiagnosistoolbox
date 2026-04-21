import copy
import sympy as sym


def create_model_copy(model_def):
    """Returns copy of user model"""
    return copy.deepcopy(model_def)


def _replace_expression(expression, substitutions):
    """Apply symbolic substitutions to SymPy expressions and leave other relation types unchanged."""
    if hasattr(expression, "xreplace"):
        return expression.xreplace(substitutions)
    return expression


def analyze_modelDef(model_def):
    """Identifies and returns replacement variables and reduced equations for the model,
    by running CSE and then undoing replacements that only depended on one original unknown variable.
    """

    # Check if input is correct type
    if not isinstance(model_def["rels"], list):
        raise TypeError("Only model_def are allowed")
    
    # Run Common Subexpression Elimination (CSE) on the model relations.
    # - replacements: list of (replacement_symbol, expression)
    # - reduced_equations: equations that use those replacement symbols
    replacements, reduced_equations = sym.cse(model_def["rels"])

    # Final list of replacements we decide to keep
    # (these are the more complex expressions worth storing as temporary variables)
    kept_replacements = []

    # Tracks each temporary symbol fully expanded using earlier replacements
    # {'symbol': expression}
    fully_expanded = {}

    # Stores temporary symbols that we want to REMOVE,
    # by directly substituting their expressions wherever they are used.
    # These are expressions that depend on only one original variable.
    # {'temporary_symbol': expression}
    replacements_to_inline = {}

    # Go through each replacement produced by CSE
    for temporary_symbol, replacement_expression in replacements:

        # Expand the expression using all previously known replacements
        # so we can see what original variables it truly depends on
        expanded_expression = _replace_expression(
            replacement_expression, fully_expanded
        )

        # Check how many variables the expression depends on
        # If it is less than 2 then we need to undo the substitution
        if len(expanded_expression.free_symbols) <= 1:

            # Mark this symbol to be removed later by substituting its value directly
            replacements_to_inline[temporary_symbol] = expanded_expression

        else:
            # This expression depends on 2 or more unknowns, keep it as a separate replacement
            # Clean it by replacing any simple symbols
            cleaned_expression = _replace_expression(
                replacement_expression, replacements_to_inline
            )

            # Store this cleaned replacement
            kept_replacements.append((temporary_symbol, cleaned_expression))

        # Save the fully expanded version for use in later steps
        fully_expanded[temporary_symbol] = expanded_expression

    # Now update the equations:
    # replace all removable temporary symbols with their actual expressions
    reduced_equations = [
        _replace_expression(equation, replacements_to_inline)
        for equation in reduced_equations
    ]

    # Return:
    # - kept_replacements: only the variables that depends on 2 or more unknown variables
    # - reduced_equations: equations where simple variables have been substituted away
    return kept_replacements, reduced_equations


def get_optimized_modelDef(input_model):
    """Returns model transformed to representation that satisfies toolbox limitations"""
    # Create model copy
    model = create_model_copy(input_model)

    replacements, reduced = analyze_modelDef(model)

    # Remove all old equations
    model["rels"] = []

    # Add all replacement variables to model
    # Add substitutions
    for replacement in replacements:
        model["x"].append(str(replacement[0]))
        model["x"].sort()

        # The replacements are pairs, e.g. (x0, x1*x2) which can be interpreted as the equation x0 = x1*x2. Since every function in rels should equal 0, x0 is moved to the other side and negated.
        model["rels"].append(-replacement[0] + replacement[1])

    #  Add replacement equations
    for element in reduced:
        model["rels"].append(element)

    return model


def logging():
    """Logs the transformations performed on the model"""
    # TODO implement function
    # använd time biblioteket för att logga tidpunkt för transformationer
    current_time = None
    pass
