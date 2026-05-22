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

    # Tracks each temporary symbol fully expanded using earlier replacements.
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


def optimize(input_model_def):
    """Transforms the model_def to a new representation with reduced expressions
    and returns a dictionary with the new model_def and a summary. The transformation
    status is printed as it progresses."""

    print("Starting model transformation")
    # Create model copy
    model = create_model_copy(input_model_def)
    print("The model was copied successfully")
    replacements, reduced = analyze_modelDef(model)
    print("Possible replacement equations have been extracted")
    # Remove all old equations
    model["rels"] = []

    # Add all replacement variables to model
    # Add substitutions
    for replacement in replacements:
        model["x"].append(str(replacement[0]))
        model["x"].sort()

        # The replacements are pairs, e.g. (x0, x1*x2) which can be interpreted as the equation x0 = x1*x2.
        # Since every function in rels should equal 0, x0 is moved to the other side and negated.
        model["rels"].append(-replacement[0] + replacement[1])

    #  Add replacement equations
    for element in reduced:
        model["rels"].append(element)

    print("Created model_def from replacement equations")
    result = {"model_def": model, "summary": create_summary(replacements, reduced)}

    return result


def get_optimized_modelDef(input_model):
    """Returns model transformed to representation that satisfies toolbox limitations"""
    return optimize(input_model)["model_def"]


def create_summary(replacements, reduced):
    """Create a pretty summary string for printing that includes replacements and reduced equations"""

    def expression_object_to_string(obj):
        """Create and return a string from expressions"""
        if isinstance(obj, str):
            return obj
        return sym.sstr(obj)

    lines = []
    lines.append("Model reduction summary")
    lines.append("=" * 80)

    lines.append("Replacements")
    lines.append("-" * 80)
    if replacements:
        for index, (left_hand_side, right_hand_side) in enumerate(
            replacements, start=1
        ):
            left_hand_side_str = expression_object_to_string(left_hand_side)
            right_hand_side_str = expression_object_to_string(right_hand_side)
            lines.append(f"{index:>3}. {left_hand_side_str} = {right_hand_side_str}")
    else:
        lines.append("  None")

    lines.append("")
    lines.append("Reduced model relations")
    lines.append("-" * 80)
    if reduced:
        for index, equation in enumerate(reduced, start=1):
            # Check if equation is derivative
            if (
                isinstance(equation, (list, tuple))
                and len(equation) == 3
                and equation[2] == "diff"
            ):
                left_hand_side_str = expression_object_to_string(equation[1])
                right_hand_side_str = expression_object_to_string(equation[0])
                lines.append(
                    f"{index:>3}. d/dt({left_hand_side_str}) = {right_hand_side_str}   [DiffConstraint]"
                )
            else:
                lines.append(f"{index:>3}. {expression_object_to_string(equation)}")
    else:
        lines.append("  None")

    lines.append("")
    lines.append("Statistics")
    lines.append("-" * 80)
    lines.append(f"  Number of replacements      : {len(replacements)}")
    lines.append(f"  Number of reduced relations : {len(reduced)}")

    return "\n".join(lines)
