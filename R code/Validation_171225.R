# ----------------------------------------------------------------------------
#
# R SCRIPT FOR VALIDATING INTERACTION ESTIMATIONS
#
# This script calculates and visualizes the percentage of "correct" (validated)
# interactions for each reconstructed network, based on expert review.
#
# Input: "Interacciones revisadas_151225.xlsx"
# Output: A bar plot showing the % of correct estimations by network.
#
# ----------------------------------------------------------------------------

# 1. Setup
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("readxl")) install.packages("readxl")
library(tidyverse)
library(readxl)

# 2. Load Data
# We read the Excel file provided by the user.
# The user mentioned "Hoja2", so we explicitly try to read that sheet.
# If "Hoja2" doesn't exist, it will fall back or error, so check your sheet names.
tryCatch({
  validation_data <- read_excel("results/Interacciones revisadas_151225.xlsx", sheet = "Hoja2")
}, error = function(e) {
  message("Sheet 'Hoja2' not found, reading the first sheet instead.")
  validation_data <- read_excel("results/Interacciones revisadas_151225.xlsx")
})

# 3. Process Data
# We need to group by network and calculate the percentage of "ok" interactions.
validation_summary <- validation_data %>%
  # Create a clean "Network" identifier combining Time and Region
  mutate(
    Network = paste0(food_web_time, "\n", food_web_region),
    # Ensure Columna1 is consistently formatted (e.g., lower case, trimmed)
    Validation = str_to_lower(str_trim(Columna1))
  ) %>%
  # Group by the unique Network ID
  group_by(Network) %>%
  # Calculate summary statistics
  summarise(
    Total_Interactions = n(),
    Correct_Estimations = sum(Validation == "ok", na.rm = TRUE),
    Percentage_Correct = (Correct_Estimations / Total_Interactions) * 100
  ) %>%
  ungroup() %>%
  # Optional: Order the networks logically if needed (e.g., by time)
  # Here we just reorder by percentage for a nicer plot, or keep alphabetical.
  # Let's order by Time (Oldest to Newest) if possible, but the names are text.
  # We can set factor levels manually for chronological order:
  mutate(
    Network = factor(Network, levels = c(
      "Pleistoceno final\nEstepa",
      "Holoceno medio\nEstepa",
      "Holoceno final\nEstepa",
      "Holoceno final\nBosque",
      "Historico\nEstepa", # Note: 'Historico' in Excel might lack accent
      "Historico\nBosque"
    ))
  )

# 4. Generate Plot
# A simple, clear bar chart showing the % correct.
p_validation <- ggplot(validation_summary, aes(x = Network, y = Percentage_Correct)) +
  geom_bar(stat = "identity", fill = "#2C3E50", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage_Correct)), 
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Validation of Reconstructed Networks",
    subtitle = "Percentage of estimated interactions confirmed as correct ('ok')",
    x = "Network (Time & Region)",
    y = "Correct Estimations (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank()
  )

# 5. Save Plot
ggsave("Network_Validation_Plot.png", p_validation, width = 8, height = 6, dpi = 300, bg = "white")

print("Validation plot generated and saved as 'Network_Validation_Plot.png'.")
print(validation_summary)
