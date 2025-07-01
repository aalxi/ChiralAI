def display_results_in_terminal(results: list[dict]):
    """
    Prints a formatted, easy-to-read summary of the feasibility analysis.
    """
    if not results:
        print("\nNo valid results to display.")
        return

    print("\n" + "="*50)
    print("🔬 CHIRALAI FEASIBILITY REPORT 🔬")
    print("="*50)

    for i, result in enumerate(results):
        name = result.get("name", "N/A")
        is_feasible = result.get("biosynthetic_feasibility", False)
        confidence = result.get("confidence", 0.0)
        host = result.get("predicted_host", "N/A")
        bottlenecks = result.get("bottlenecks", [])

        status_emoji = "✅ Feasible" if is_feasible else "❌ Not Feasible"
        
        print(f"\n--- Result {i+1}: {name} ---")
        print(f"  - Status: {status_emoji}")
        print(f"  - Confidence: {confidence:.0%}")
        print(f"  - Host: {host}")
        
        # Only show bottlenecks if they exist
        if bottlenecks and bottlenecks[0] != "None identified":
            print("  - Bottlenecks:")
            for bn in bottlenecks:
                print(f"    • {bn}")

    print("\n" + "="*50)
    print("Full details saved in the 'output' directory.")